import pandas as pd 


#colnames = ['chr','pos','dbsnp','MAF','ref','alt','A','C','G','T','RD','Ars','Crs','Grs','Trs','af','cov']
colnames = ['chr','pos']

df1 = pd.read_csv('prova1', names=colnames, sep='\t')
df2 = pd.read_csv('prova2', names=colnames, sep='\t')
df3 = pd.read_csv('prova3', names=colnames, sep='\t')

print df1
set1 = df1.apply(set, axis =1)
print set1

Afilt=[df1, df2, df3]
'''
dfA_merge = Afilt[0]
for df in Afilt[1:]:       
    dfA_merge = pd.merge(dfA_merge, df, on='dbsnp', how='inner')

#print dfA_merge

dfA_merge.to_csv('mergeprova.csv', sep='\t')

#print pd.concat(Afilt).drop_duplicates().reset_index(drop=True)
'''