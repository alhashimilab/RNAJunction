#!/usr/bin/python

import json
import os
import sys
import time
import pandas as pd
import learnna_json as lna_json
from pdblib.base import *
from commontool import read, readchar

def addTER(filename):  #Add TER in the dssr-bulge file 
    strings = read(filename)  #inport bulge file into string list
    chars = readchar(filename)  #import bulge file into char list
    ters = []  #line index that should add ter
    count = 0  #initialize count of residue id
    for i in range(len(strings)):
        if strings[i][0] == 'REMARK' and strings[i+1][0] == 'REMARK' and strings[i+2][0] == 'REMARK':
            count = 1;
        if strings[i][0] in ['ATOM','HETATM'] and strings[i+1][0] in ['ATOM','HETATM']:
            if chars[i][21] != chars[i+1][21] and count == 1:  #read column of resi
                count = 0  #count + 1 if the resi change
                ters.append(i)
    for i in range(len(ters)):
        os.system("sed -i '%d i\TER' %s"%(ters[i]+i+2,filename))
    return ters

inpdir = './RawBulge'
jsondir = './Json'
oupdir = './Bulge'
curdir = os.getcwd()
bulgedic = {}  #dictionary to store pdb -- bulge information
pdblist = []  #list to store pdb file names

#Walk pdb_directory and get list of pdb names
for file in read('Pdbinfo/All_nmr.txt'):
    pdblist.append(file[0])

pdblist.sort()  #sort the pdblist

tic = time.time()  #timer: start

if len(os.listdir(oupdir)) != 0:
    print ("Files already in output folder! Please check.")
    sys.exit()

for idx, pdbid in enumerate(pdblist):
    print ("--- Working on [%s] (%d of %d) ---"%(pdbid,idx+1,len(pdblist)))
    pdb = pdbid + ".pdb"
    bulgedic[pdbid] = ''
    #Step1: Check the dssr-bulges file in RawBulge folder
    ipdbf = os.path.join(inpdir,pdb)
    jsonf = os.path.join(jsondir,pdb.replace(".pdb",".json"))
    if pdb not in os.listdir(inpdir):  #check if we already process the pdbfile
        continue
    #Step2: Read the json file
    najson = lna_json.NA_JSON()  #initialize class objects
    with open(jsonf) as json_data:  #read each json file
        data = json.load(json_data)
    najson.set_json(data)  #pass json file to class pbject
    najson.read_idx()  #set index from own json file
    bulgesinfo = najson.json_file['bulges']
    #Step4: Add TER between two strands for each bulge in bulge file
    ter = addTER(ipdbf)
    #Step5: Split bulge file into each bulge using pdblib
    bulges = Pdb(ipdbf)
    for i, md in enumerate(bulges.mds):
        if len(md.segs) > 2:
            print(">>>Abort! There are more than 2 chains: %s [%d/%d]"%(pdb,idx+1,len(pdblist)))
            sys.exit()
        bulgeconf = bulgesinfo[i]['bridging_nts']
        bulgenum = max(bulgeconf[0],bulgeconf[1])
        opdbf = os.path.join(oupdir,pdbid+"_"+str(i)+"_"+str(bulgenum)+".pdb")  #pdbid-bulgeid-bulgenum.pdb
        bulgedic[pdbid] += " "+pdbid+"_"+str(i)+"_"+str(bulgenum)
        md.write("%s"%opdbf)
    print(bulgedic[pdbid])

toc = time.time()  #timer: end
print(">>>>>> Total time used %f >>>>>>"%(toc-tic))

bulge_df = pd.DataFrame()
bulge_df['pdb'] = bulgedic.keys()
bulge_df['bulge'] = bulgedic.values()
bulge_df = bulge_df.sort(['pdb']).reset_index(drop=True)
bulge_df.to_csv("Bulge_nmr.csv",index=False)
