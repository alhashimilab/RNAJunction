#!/usr/bin/python

import pandas as pd

df00 = pd.read_csv('PDB_match_conditions.csv')
df01 = pd.read_csv('PDB_match_nucleosome.csv')
df02 = pd.read_csv('PDB_match_ligand.csv',keep_default_na=False, na_values=[''])
df11 = pd.read_csv('Stem_crystal.csv')
df12 = pd.read_csv('Bulge_crystal.csv')
df13 = pd.read_csv('Iloop_crystal.csv')
df14 = pd.read_csv('Hairpin_crystal.csv')


df = pd.concat([df11.pdb, df00.molecule_type, df00.classification, df00.resolution, df01.Name, df02.Ligand, df02.Metal, df02.Solvent, df11.stem, df12.bulge, df13.iloop, df14.hairpin],axis=1)

df.to_csv('Final_crystal.csv',index=False,header = ['pdb','mmtype','mmclass','reso','name','ligand','metal','solvent','stem','bulge','iloop','hairpin'],keep_default_na=False, na_values=[''])

