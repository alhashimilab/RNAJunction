#!/usr/bin/python

import json
import os
import sys
import time
import pandas as pd
import numpy as np
import learnna_json as lna_json
from pdblib.abg import *
from commontool import read, readchar


oupname = 'IloopTable_Crystal_2bp_dawn'
inpdir = './Iloop'
jsondir = './Json'
jsonlist = []
df = pd.read_csv('./Pdbinfo/Final_crystal.csv')
ilooppdb = df[df.iloop.notnull()]
pdblist = ilooppdb.pdb.tolist()
ilooplist = [['pdbid','mmtype','mmclass','reso','iloopid','iloopnum1','iloopnum2','chid1','chid2','seq_short1','seq_short2','seq_long1','seq_long2','modify','abg_seq1','abg_seq2','dotbra_1','dotbra_2','RMSD1','RMSD2','alpha_h','beta_h','gamma_h']]
helixlen = 2
biglist = [a[0] for a in read('biglist.txt')]
fnhlx = './iAformRNA.pdb'
hlx = Mol(fnhlx)

tic = time.time()  #start time
for idx, pdbid in enumerate(pdblist):
    #print ("--- Working on [%s] (%d of %d) ---"%(pdbid,idx+1,len(pdblist)))
    #Step1: Read the json file
    abg_seq1_set = set()
    abg_seq2_set = set()
    dotbra_1_set = set()
    dotbra_2_set = set()
    mmtype = ilooppdb.mmtype[ilooppdb.pdb == pdbid].values[0]
    mmclass = ilooppdb.mmclass[ilooppdb.pdb == pdbid].values[0]
    reso = ilooppdb.reso[ilooppdb.pdb == pdbid].values[0]
    jsonf = os.path.join(jsondir,pdbid+'.json') 
    najson = lna_json.NA_JSON()  #initialize class objects
    with open(jsonf) as json_data:  #read each json file
        data = json.load(json_data)
    najson.set_json(data)  #pass json file to class pbject
    najson.read_idx()  #set index from own json file
    #Step2: Store the index of each nucleotide into a dictionary
    nts_idx = {}  #nts index dictionary with nt_id as key
    for i, nt in enumerate(najson.json_file['nts']):
        nts_idx[nt['nt_id']] = i
    if len(nts_idx) != len(najson.json_file['nts']):  #error: there are nucleotides labeled same
        sys.exit()
    nts = najson.json_file['nts']
    iloops = najson.json_file['iloops']
    #Step3: loop all the iloop in one json file
    for i, iloop in enumerate(iloops):
        bridging_nts = iloop['bridging_nts']
        iloopnum1 = bridging_nts[0]
        iloopnum2 = bridging_nts[1]
        iloopid = pdbid + '_' + str(i) + '_' + str(iloopnum1) + '_' + str(iloopnum2)
        iloopf = os.path.join(inpdir,iloopid+'.pdb')
        #check if pdb file exist 
        if os.path.isfile(iloopf) == False:
            continue
        #get segment from iloop pdb file
        iloopmol = Mol(iloopf)
        if len(iloopmol.segs) == 2:
            iloopseg1 = iloopmol.segs[0]
            iloopseg2 = iloopmol.segs[1]
        else:
            print("ERROR: %10s Incorrect segment number"%iloopid)
            sys.exit()
        bridging_nts = iloop['bridging_nts']
        bridges = iloop['bridges']
        len_nts_1 = iloopnum1 + 2
        len_nts_2 = iloopnum2 + 2
        seq_short1 = bridges[0]['nts_short']
        seq_short2 = bridges[1]['nts_short']
        seq_long1 = bridges[0]['nts_long']
        seq_long2 = bridges[1]['nts_long']
        nts_long = iloop['nts_long'].split(',')
        nts_long_1 = nts_long[0:len_nts_1]
        nts_long_2 = nts_long[len_nts_1:]
        #get edge nts_id in json file
        nts_1_0 = nts_long_1[0]
        nts_1_1 = nts_long_1[-1]
        nts_2_0 = nts_long_2[0]
        nts_2_1 = nts_long_2[-1]
        #print "nts: " + nts_1_0 + " " + nts_1_1 + " " + nts_2_0 + " " + nts_2_1
        #assign chain id
        chid1 = iloopseg1.chid
        chid2 = iloopseg2.chid
        chid_1_0 = nts[nts_idx[nts_1_0]]['chain_name']
        chid_1_1 = nts[nts_idx[nts_1_1]]['chain_name']
        chid_2_0 = nts[nts_idx[nts_2_0]]['chain_name']
        chid_2_1 = nts[nts_idx[nts_2_1]]['chain_name']
#        if (chid1 != chid_1_0) or (chid1 != chid_1_1) or (chid2 != chid_2_0) or (chid2 != chid_2_1):
#            print("ERROR: chain id is inconsistent!")
#            sys.exit()
        #assign residue name
        resn_1_0 = nts[nts_idx[nts_1_0]]['nt_name']
        resn_1_1 = nts[nts_idx[nts_1_1]]['nt_name']
        resn_2_0 = nts[nts_idx[nts_2_0]]['nt_name']
        resn_2_1 = nts[nts_idx[nts_2_1]]['nt_name']
        #assign residue id & icode
        resic_1_0 = nts[nts_idx[nts_1_0]]['nt_id'].split('.')[1].replace(resn_1_0,"").replace("/","")
        resic_1_1 = nts[nts_idx[nts_1_1]]['nt_id'].split('.')[1].replace(resn_1_1,"").replace("/","")
        resic_2_0 = nts[nts_idx[nts_2_0]]['nt_id'].split('.')[1].replace(resn_2_0,"").replace("/","")
        resic_2_1 = nts[nts_idx[nts_2_1]]['nt_id'].split('.')[1].replace(resn_2_1,"").replace("/","")
        #assign icode: check the existence of icode in resic, if exist, icode = icode; else, icode = " "
        if resic_1_0.find("^") != -1:
            resi_1_0 = int(resic_1_0.split("^")[0])
            icode_1_0 = str(resic_1_0.split("^")[1])
        else:
            resi_1_0 = int(resic_1_0)
            icode_1_0 = " "
        if resic_1_1.find("^") != -1:
            resi_1_1 = int(resic_1_1.split("^")[0])
            icode_1_1 = str(resic_1_1.split("^")[1])
        else:
            resi_1_1 = int(resic_1_1)
            icode_1_1 = " "
        if resic_2_0.find("^") != -1:
            resi_2_0 = int(resic_2_0.split("^")[0])
            icode_2_0 = str(resic_2_0.split("^")[1])
        else:
            resi_2_0 = int(resic_2_0)
            icode_2_0 = " "
        if resic_2_1.find("^") != -1:
            resi_2_1 = int(resic_2_1.split("^")[0])
            icode_2_1 = str(resic_2_1.split("^")[1])
        else:
            resi_2_1 = int(resic_2_1)
            icode_2_1 = " "
#        if icode_1_0 != " " or icode_1_1 != " " or icode_2_0 != " " or icode_2_1 != " ":
#            print("WARNING: icode exists in edge residue, needs your attention!")
        #print "resi: " + str(resi_1_0) + " " + str(resi_1_1) + " " + str(resi_2_0) + " " + str(resi_2_1)
        #print "icode: " + icode_1_0 + " " + icode_1_1 + " " + icode_2_0 + " " + icode_2_1
        res_1_0 = getres(iloopseg1,resi_1_0,icode_1_0)
        res_1_1 = getres(iloopseg1,resi_1_1,icode_1_1)
        res_2_0 = getres(iloopseg2,resi_2_0,icode_2_0)
        res_2_1 = getres(iloopseg2,resi_2_1,icode_2_1)
        reses_1 = getreses(iloopseg1)
        reses_2 = getreses(iloopseg2)
        index_1_0 = reses_1.index(res_1_0)
        index_1_1 = reses_1.index(res_1_1)
        index_2_0 = reses_2.index(res_2_0)
        index_2_1 = reses_2.index(res_2_1)
        # print "index: " + str(index_1_0) + " " + str(index_1_1) + " " + str(index_2_0) + " " + str(index_2_1)
        # check the rationality of index
        if index_1_0 >= index_1_1 or index_2_0 >= index_2_1:
            print index_1_0
            print index_1_1
            print index_2_0
            print index_2_1
            print("ERROR: %10s residue id is not read properly using pdblib"%iloopid)
            sys.exit()
        if len(reses_1) == len(reses_2) and iloopnum1 != iloopnum2:
            print("ERROR: %10s same residue number on both strand"%iloopid)
            sys.exit()
        # check boundary of res index in reses list
        if index_1_0-(helixlen-1) < 0 or index_1_1+(helixlen-1) >= len(reses_1) or index_2_0-(helixlen-1) < 0 or index_2_1+(helixlen-1) >= len(reses_2):
            print("Skip!: %10s index_1_0: %d index_1_1: %d ilooplen1: %d index_2_0: %d index_2_1: %d ilooplen2: %d"%(iloopid,index_1_0,index_1_1,len(reses_1),index_2_0,index_2_1,len(reses_2)))
            continue
        res_long_1 = reses_1[index_1_0-(helixlen-1):index_1_1+(helixlen-1)+1]
        res_long_2 = reses_2[index_2_0-(helixlen-1):index_2_1+(helixlen-1)+1]
        # assign h1 and h2 for abg calculation
        if (len(nts_long_1) >= len(nts_long_2)):
            mresl1 = res_long_1[0:helixlen] + res_long_2[-(helixlen):]
            mresl2 = res_long_1[-(helixlen):] + res_long_2[0:helixlen]
        elif (len(nts_long_1) < len(nts_long_2)):
            mresl1 = res_long_2[0:helixlen] + res_long_1[-(helixlen):]
            mresl2 = res_long_2[-(helixlen):] + res_long_1[0:helixlen]
        # find modification
        modify = []
        # counting icode number in the res_long list
        icode_num_1 = 0
        icode_num_2 = 0
        for res in res_long_1: 
            if res.name not in ['DA','DG','DC','DT','A','G','C','U']:
                modify_nt = res.chid+'.'+res.name+'/'+str(res.resi)+res.icode
                modify.append(modify_nt.replace(" ",""))
            if res.icode != " ":
                icode_num_1 = icode_num_1 + 1
        for res in res_long_2:
            if res.name not in ['DA','DG','DC','DT','A','G','C','U']:
                modify_nt = res.chid+'.'+res.name+'/'+str(res.resi)+res.icode
                modify.append(modify_nt.replace(" ",""))
            if res.icode != " ":
                icode_num_2 = icode_num_2 + 1
        if icode_num_1 != 0:
            print("Icode: %10s icodenum1: %4d resi0: %4d resi-1: %4d res_long_len: %4d "%(iloopid,icode_num_1,res_long_1[0].resi,res_long_1[-1].resi,len(res_long_1)))
        if icode_num_2 != 0:
            print("Icode: %10s icodenum2: %4d resi0: %4d resi-1: %4d res_long_len: %4d "%(iloopid,icode_num_2,res_long_2[0].resi,res_long_2[-1].resi,len(res_long_2)))
        modify = ','.join(modify)
        # check if residue id is consistent with the res_long length considering icode
        #if res_long_1[0].resi + len(res_long_1) - icode_num_1 - 1 != res_long_1[-1].resi or res_long_2[0].resi + len(res_long_2) - icode_num_2 - 1 != res_long_2[-1].resi:
        #   print("ERROR: %10s The residue id is not consistent considering icode!"%iloopid)
        #   sys.exit()
        abg_seq1 = chid1 + "." + str(res_long_1[0].resi) + str(res_long_1[0].icode) + ":" + str(res_long_1[-1].resi) + str(res_long_1[-1].icode)
        abg_seq2 = chid2 + "." + str(res_long_2[0].resi) + str(res_long_2[0].icode) + ":" + str(res_long_2[-1].resi) + str(res_long_2[-1].icode)
        abg_seq1 = abg_seq1.replace(" ","")
        abg_seq2 = abg_seq2.replace(" ","")
        dotbra_1 = '(((('
        dotbra_2 = '))))'
        dotbra_1 = "'" + dotbra_1[0:(helixlen)] + "."*(len_nts_1 - 2) + dotbra_1[(helixlen):] + "'"
        dotbra_2 = "'" + dotbra_2[0:(helixlen)] + "."*(len_nts_2 - 2) + dotbra_2[(helixlen):] + "'"
        # deal with the redundancy with pdb in the biglist
        if pdbid in biglist:
            if abg_seq1 in abg_seq1_set and abg_seq2 in abg_seq2_set and dotbra_1 in dotbra_1_set and dotbra_2 in dotbra_2_set: 
                continue
            else:
                abg_seq1_set.add(abg_seq1)
                abg_seq2_set.add(abg_seq2)
                dotbra_1_set.add(dotbra_1)
                dotbra_2_set.add(dotbra_2)
                rmsd1, rmsd2, alpha_h, beta_h, gamma_h = getabgA2(hlx,mresl1,mresl2)
                ilooplist.append([pdbid,mmtype,mmclass,reso,iloopid,iloopnum1,iloopnum2,chid1,chid2,seq_short1,seq_short2,seq_long1,seq_long2,modify,abg_seq1,abg_seq2,dotbra_1,dotbra_2,rmsd1,rmsd2,-alpha_h,beta_h,-gamma_h])
        else:    
            rmsd1, rmsd2, alpha_h, beta_h, gamma_h = getabgA2(hlx,mresl1,mresl2)
            ilooplist.append([pdbid,mmtype,mmclass,reso,iloopid,iloopnum1,iloopnum2,chid1,chid2,seq_short1,seq_short2,seq_long1,seq_long2,modify,abg_seq1,abg_seq2,dotbra_1,dotbra_2,rmsd1,rmsd2,-alpha_h,beta_h,-gamma_h])
        print iloopid + " " + str(nts_long_1) + " " + str(nts_long_2)
        

toc = time.time()  #end time
print(">>>>>> Total time used %f >>>>>>"%(toc-tic))

iloop_df = pd.DataFrame(ilooplist[1:], columns=ilooplist[0])
#iloop_df = iloop_df.loc[(iloop_df.mmtype == 'RNA') | (iloop_df.mmtype == 'Protein#RNA')]
iloop_df.to_csv("%s.csv"%oupname,index=False)
