import pandas as pd
from commontool import read

ns_df = pd.read_csv('Nucleosome_query.csv')

nsl = ns_df['PDB ID'].tolist()
print nsl

pdblist = [a[0].upper() for a in read('All_crystal.txt')]
ns_dic = dict([(key, None) for key in pdblist])

for ns in nsl:
    ns_dic[ns] = 'Nucleosome'

newdf = pd.DataFrame(ns_dic.items(), columns=['PDBID','Name'])
newdf = newdf.sort(['PDBID'])
newdf.to_csv('PDB_match_nucleosome.csv',index=False)



