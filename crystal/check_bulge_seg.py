#!/usr/bin/python

import json
import os
import sys
import time
import pandas as pd
import numpy as np
import learnna_json as lna_json
from pdblib.base import *
from commontool import read, readchar


inpdir = './Bulge'
pdblist = os.listdir(inpdir)
pdblist.sort()

tic = time.time()  #start time
for idx, pdbid in enumerate(pdblist):
    #print ("--- Working on [%s] (%d of %d) ---"%(pdbid,idx+1,len(pdblist)))
    #Step1: Read the json file
    bulgef = inpdir + '/' + pdbid
    bulgemol = Mol(bulgef)
    if len(bulgemol.segs) == 2:
        for i in range(len(bulgemol.segs[0].reses)-1):
            if bulgemol.segs[0].reses[i].resi + 1 != bulgemol.segs[0].reses[i+1].resi:
                print "Abnormal residue index: " + pdbid
                break
        for i in range(len(bulgemol.segs[1].reses)-1):
            if bulgemol.segs[1].reses[i].resi + 1 != bulgemol.segs[1].reses[i+1].resi:
                print "Abnormal residue index: " + pdbid
                break
    elif len(bulgemol.segs) > 2:
        print "ERROR: More than 2 segment: " + pdbid
        sys.exit()
    elif len(bulgemol.segs) == 1:
        print "Only 1 segment: " + pdbid
        newmol = Mol()
        newmol.segs = [Segment(),Segment()]
        oldreses = bulgemol.segs[0].reses
        ter = -1
        for i in range(len(oldreses)-1):
            if oldreses[i].resi + 1 != oldreses[i+1].resi and oldreses[i].icode == ' ' and oldreses[i+1].icode == ' ':
                ter = i
        if ter == -1:
            continue
        else:
            newmol.segs[0].reses = oldreses[0:ter+1]
            newmol.segs[1].reses = oldreses[ter+1:]
            newmol.segs[0].chid = bulgemol.segs[0].chid
            newmol.segs[1].chid = bulgemol.segs[0].chid
            newmol.write(bulgef)
