#!/usr/bin/python

import json
import os
import sys
import time
import pandas as pd
import learnna_json as lna_json
from pdblib.base import *
from commontool import read, readchar

def addTER(filename):  #Add TER in the dssr-iloop file 
    strings = read(filename)  #inport iloop file into string list
    chars = readchar(filename)  #import iloop file into char list
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

inpdir = './RawIloop'
jsondir = './Json'
oupdir = './Iloop'
curdir = os.getcwd()
iloopdic = {}  #dictionary to store pdb -- iloop information
pdblist = []  #list to store pdb file names

#Walk pdb_directory and get list of pdb names
for file in read('Pdbinfo/All_crystal.txt'):
    pdblist.append(file[0])

pdblist.sort()  #sort the pdblist

tic = time.time()  #timer: start

if len(os.listdir(oupdir)) != 0:
    print ("Files already in output folder! Please check.")
    sys.exit()

for idx, pdbid in enumerate(pdblist):
    print ("--- Working on [%s] (%d of %d) ---"%(pdbid,idx+1,len(pdblist)))
    pdb = pdbid + ".pdb"
    iloopdic[pdbid] = ''
    #Step1: Check  the dssr-iloops file in RawIloop folder
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
    iloopsinfo = najson.json_file['iloops']
    #Step3: Add TER between two strands for each iloop in iloop file
    ter = addTER(ipdbf)
    #Step4: Split iloop file into each iloop using pdblib
    iloops = Pdb(ipdbf)
    for i, md in enumerate(iloops.mds):
        if len(md.segs) > 2:
            print(">>>Abort! There are more than 2 chains: %s [%d/%d]"%(pdb,idx+1,len(pdblist)))
            sys.exit()
        iloopconf = iloopsinfo[i]['bridging_nts']
        iloopnum1 = iloopconf[0]
        iloopnum2 = iloopconf[1]
        opdbf = os.path.join(oupdir,pdbid+"_"+str(i)+"_"+str(iloopnum1)+"_"+str(iloopnum2)+".pdb")  #pdbid-iloopid-iloopnum.pdb
        iloopdic[pdbid] += " "+pdbid+"_"+str(i)+"_"+str(iloopnum1)+"_"+str(iloopnum2)
        md.write("%s"%opdbf)
    print(iloopdic[pdbid])

toc = time.time()  #timer: end
print(">>>>>> Total time used %f >>>>>>"%(toc-tic))

iloop_df = pd.DataFrame()
iloop_df['pdb'] = iloopdic.keys()
iloop_df['iloop'] = iloopdic.values()
iloop_df = iloop_df.sort(['pdb']).reset_index(drop=True)
iloop_df.to_csv("Iloop_crystal.csv",index=False)
