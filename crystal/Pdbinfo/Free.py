import pandas as pd

df = pd.read_csv('Final_crystal.csv',keep_default_na=False, na_values=[''])

newdf = df.loc[(df.mmtype == 'RNA') & ~((df.ligand.isnull()) & (df.metal.isnull()) & (df.solvent.isnull()))]

newdf.to_csv('RNA_crystal.csv')
