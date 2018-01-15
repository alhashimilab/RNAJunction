import pandas as pd
from commontool import read

modify = [m[0] for m in read('modify.txt')]
ligand = [m[0] for m in read('ligand.txt')]
solvent = [m[0] for m in read('solvent.txt')]
metal = [m[0] for m in read('metal.txt')]
protein = [m[0] for m in read('protein.txt')]

df = pd.read_csv('PDB_match_conditions.csv',keep_default_na=False, na_values=['_'])

for i in range(len(df)):
    pdb_id = df.ix[i]['pdb_id']
    ll = df.ix[i]['ligandid'].split(" | ")
    for l in ll:
        if l not in modify and l not in ligand and l not in solvent and l not in metal and l not in protein and l != 'null' and l != 'UNX':
            print pdb_id + " " + l

