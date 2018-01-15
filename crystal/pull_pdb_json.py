#!/usr/bin/python

import os
import sys
import time
import pandas as pd
from pdblib.base import *
from commontool import read, readchar

inpdir = './Crystal'
jsondir = './Json'
oupdir_helice = './RawHelice' # raw stem fragment directory
oupdir_stem = './RawStem' # raw stem fragment directory
oupdir_hairpin = './RawHairpin' # raw hairpin fragment directory
oupdir_bulge = './RawBulge' # raw bulge fragement directory
oupdir_iloop = './RawIloop' # raw iloop fragment directory
curdir = os.getcwd()
pdblist = []  # list to store pdb file names
jsonlist = []  # list to store existing json file names

# loop pdb_directory and get list of pdb names
for file in read('./Pdbinfo/All_crystal.txt'):
    pdblist.append(file[0])

pdblist.sort()  #sort the pdblist

# loop json_directory to check if the pdb exist in the json file
for file in os.listdir(jsondir):
    jsonlist.append(file)

jsonlist.sort()  #sort the pdblist

tic = time.time()  #timer: start

for idx, pdbid in enumerate(pdblist):
    print ("--- Working on [%s] (%d of %d) ---"%(pdbid,idx+1,len(pdblist)))
    pdb = pdbid + ".pdb"
    #Step1: Create output json file and dssr files
    ipdbf = os.path.join(inpdir,pdb)
    jsonf = pdb.replace(".pdb",".json")
    if jsonf in jsonlist:  #check if we already process the pdbfile
        continue
    os.system("x3dna-dssr --more --loop=with-stems --json --symm -i=%s -o=%s"%(ipdbf,jsonf))
    #Step2: Check if there is output dssr-helice file
    if "dssr-helices.pdb" in os.listdir(curdir):
        os.system("cp dssr-helices.pdb %s"%os.path.join(oupdir_helice,pdb))
    ##Step3: Check if there is output dssr-stem file
    if "dssr-stems.pdb" in os.listdir(curdir):
        os.system("cp dssr-stems.pdb %s"%os.path.join(oupdir_stem,pdb))
    ##Step4: Check if there is output dssr-hairpin file
    if "dssr-hairpins.pdb" in os.listdir(curdir):
        os.system("cp dssr-hairpins.pdb %s"%os.path.join(oupdir_hairpin,pdb))
    ##Step5: Check if there is output dssr-bulge file
    if "dssr-bulges.pdb" in os.listdir(curdir):
        os.system("cp dssr-bulges.pdb %s"%os.path.join(oupdir_bulge,pdb))
    ##Step6: Check if there is output dssr-iloop file
    if "dssr-iloops.pdb" in os.listdir(curdir):
        os.system("cp dssr-iloops.pdb %s"%os.path.join(oupdir_iloop,pdb))
    os.system("cp %s %s"%(jsonf,os.path.join(jsondir,jsonf)))
    os.system("x3dna-dssr --cleanup")
    os.system("rm %s"%jsonf)

toc = time.time()  #timer: end
print(">>>>>> Total time used %f >>>>>>"%(toc-tic))

