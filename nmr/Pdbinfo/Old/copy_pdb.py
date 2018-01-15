import pandas as pd
import os

df = pd.read_csv('RNA_nmr.csv')
pdblist = df.pdb.tolist()

for pdb in pdblist:
    os.system('sudo cp /mnt/hs189/NAfinder_2.0_nmr/Json/%s.json /home/hs189/Others/yuze/HbondFinder/Json_nmr/'%pdb)
    #os.system('sudo cp /mnt/hs189/NAfinder_2.0_nmr/Nmr/%s.pdb /home/hs189/Others/yuze/PDB_class/Nmr/'%pdb)
