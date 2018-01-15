import time
import pandas as pd
import numpy as np
from commontool import read

df = pd.read_csv('./Ligand_query.csv',keep_default_na=False, na_values=[''])
pdblist = [a[0].upper() for a in read('All_crystal.txt')]
metallist = [a[0] for a in read('metal.txt')]
solventlist = [a[0] for a in read('solvent.txt')]

ligand_dic = dict([(key, None) for key in pdblist])
metal_dic = dict([(key, None) for key in pdblist])
solvent_dic = dict([(key, None) for key in pdblist])
modify_dic = dict([(key, None) for key in pdblist])

tic = time.time()  #start time

for i in range(len(df)):
    ligand_id = df.ix[i]['Ligand ID']
    print ("--- Working on [%s] (%d of %d) ---"%(ligand_id,i+1,len(df)))
    #assign ligand
    if ligand_id in metallist:
        dic = metal_dic
    elif ligand_id in solventlist:
        dic = solvent_dic
    else:
        dic = ligand_dic
    pdb_ligand = df.ix[i]['Instance PDB IDs as Free Ligand'].split(',')
    for pl in (pl for pl in pdb_ligand if pl != ""):
        if pl in pdblist:
            if dic[pl] == None:
                dic[pl] = []
            if ligand_id not in dic[pl]:
                dic[pl].append(ligand_id)
    #assign modify
    pdb_modify = df.ix[i]['Instance PDB IDs as part of a polymer'].split(',')
    for ml in (ml for ml in pdb_modify if ml != ""):
        if ml in pdblist:
            if modify_dic[ml] == None:
                modify_dic[ml] = []
            if ligand_id not in modify_dic[ml]:
                modify_dic[ml].append(ligand_id)

for key in pdblist:
    if ligand_dic[key] != None:
        ligand_dic[key] = " ".join(str(x) for x in ligand_dic[key])
    if metal_dic[key] != None:
        metal_dic[key] = " ".join(str(x) for x in metal_dic[key])
    if solvent_dic[key] != None:
        solvent_dic[key] = " ".join(str(x) for x in solvent_dic[key])
    if modify_dic[key] != None:
        modify_dic[key] = " ".join(str(x) for x in modify_dic[key])

toc = time.time()  #timer: end
print(">>>>>> Total time used %f >>>>>>"%(toc-tic))

df_ligand = pd.DataFrame(ligand_dic.items(), columns=['PDBID','Ligand'])
df_metal = pd.DataFrame(metal_dic.items(), columns=['PDBID','Metal'])
df_solvent = pd.DataFrame(solvent_dic.items(), columns=['PDBID','Solvent'])
df_modify = pd.DataFrame(modify_dic.items(), columns=['PDBID','Modify'])

df_ligand = df_ligand.sort(['PDBID'])
df_metal = df_metal.sort(['PDBID'])
df_solvent = df_solvent.sort(['PDBID'])
df_modify = df_modify.sort(['PDBID'])

newdf = pd.concat([df_ligand.PDBID, df_ligand.Ligand, df_metal.Metal, df_solvent.Solvent, df_modify.Modify],axis=1)

newdf.to_csv('PDB_match_ligand.csv',index=False)


