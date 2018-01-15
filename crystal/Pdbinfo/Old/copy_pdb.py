import pandas as pd
import os

df = pd.read_csv('RNA_crystal.csv')
pdblist = df.pdb.tolist()

for pdb in pdblist:
    os.system('sudo cp /mnt/hs189/NAfinder_2.0_crystal/Json/%s.json /home/hs189/Others/yuze/HbondFinder/Json_crystal/'%pdb)
    #os.system('sudo cp /mnt/hs189/NAfinder_2.0_crystal/Crystal/%s.pdb /home/hs189/Others/yuze/PDB_class/Crystal/'%pdb)
