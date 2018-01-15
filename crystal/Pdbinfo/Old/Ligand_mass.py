import pandas as pd

df = pd.read_csv('Ligand_query.csv',keep_default_na=False, na_values=[''])

print df.columns

newdf = pd.DataFrame()
newdf['Ligand ID'] = df['Ligand ID']
newdf['Ligand MW'] = df['Ligand MW']
newdf['Ligand Formula'] = df['Ligand Formula']
newdf['Ligand Name'] = df['Ligand Name']
newdf['Ligand SMILES'] = df['Ligand SMILES']

print newdf
newdf.to_csv('Ligand_mass_crystal.csv',index=False)
