import pandas as pd 
import glob
import os 
import plotly.graph_objects as go
from itertools import cycle
from supervenn import supervenn
import matplotlib.pyplot as plt 
import re



colnames = ['chr','pos','dbsnp','MAF','ref','alt','A','C','G','T','RD','Ars','Crs','Grs','Trs','af','cov']
'''
#colnames = ['chr','pos']


df1 = pd.read_csv('prova1', names=colnames, sep='\t')
df2 = pd.read_csv('prova2', names=colnames, sep='\t')
#df3 = pd.read_csv('prova3', names=colnames, sep='\t')

#set1 = df1.apply(set, axis =1)
#print set1

#Afilt=[df1, df2, df3]

#dfA_merge = Afilt[0]
#for df in Afilt[1:]:       
 #   dfA_merge = pd.merge(dfA_merge, df, on='dbsnp', how='inner')

#print dfA_merge

#dfA_merge.to_csv('mergeprova.csv', sep='\t')

#print pd.concat(Afilt).drop_duplicates().reset_index(drop=True)

intersect=pd.merge(df1, df2, on="dbsnp", how='inner')

print intersect

opt=((intersect['af_x'] >= 0.1) & (intersect['af_x'] <= 0.9)) | ((intersect['af_y'] >= 0.1) & (intersect['af_y'] <= 0.9))

rslt_df=intersect[opt]

print rslt_df
'''
mocks_n_cov10=[]

for file in glob.glob('cov_filtering/notre/mock_TOT*.10filtered'):
	name = file.split("/")
	supp = name[2].split(".")
	df = 'df_' + supp[0]
	print df 
	#df = pd.read_csv(file, names=colnames, sep='\t')

	#print(df)
	#mocks_n_cov10.append(df)
