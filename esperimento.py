import pandas as pd 
import glob
import os 
import plotly.graph_objects as go
from itertools import cycle
from supervenn import supervenn
import matplotlib.pyplot as plt 
import re
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

pd.options.mode.chained_assignment = None  # default='warn'
colnames = ['chr','pos','dbsnp','MAF','ref','alt','A','C','G','T','RD','Ars','Crs','Grs','Trs','af','cov']
'''
#colnames = ['chr','pos']


df1 = pd.read_csv('prova1', sep='\t')
df2 = pd.read_csv('prova2', sep='\t')
#df3 = pd.read_csv('prova3', sep='\t')

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

mocks_n_cov10=[]

for file in glob.glob('cov_filtering/notre/mock_TOT*.10filtered'):
	name = file.split("/")
	supp = name[2].split(".")
	df = 'df_' + supp[0]
	print df 
	#df = pd.read_csv(file, sep='\t')

	#print(df)
	#mocks_n_cov10.append(df)

a = [1,2,3]
arr=[]
for x in range(3):
	arr.append(a[x % len(a)])
	arr.append(a[(x+1) % len(a)])
	print arr
	arr=[]



mock_TOT_1='cov_filtering/notre/mock_TOT_1_S15.final.PILEUP.ASEQ.10filtered'
mock_TOT_3='cov_filtering/notre/mock_TOT_3_S14.final.PILEUP.ASEQ.10filtered'
mock_TOT_6='cov_filtering/notre/mock_TOT_6_S13.final.PILEUP.ASEQ.10filtered'


mock_TOT_1='af_filtering/local1/mock_TOT_1_S15.final.PILEUP.ASEQ.Dfiltered'
mock_TOT_3='af_filtering/local1/mock_TOT_3_S14.final.PILEUP.ASEQ.Dfiltered'
mock_TOT_6='af_filtering/local1/mock_TOT_6_S13.final.PILEUP.ASEQ.Dfiltered'

#creation of dataframes
df1TR = pd.read_csv(mock_TOT_1, sep='\t')
df3TR = pd.read_csv(mock_TOT_3, sep='\t')
df6TR = pd.read_csv(mock_TOT_6, sep='\t')

inter = pd.merge(df1TR, df3TR, on="dbsnp", how='inner')

df=inter[((inter['af_x'] >= 0.1) & (inter['af_x'] <= 0.9)) | ((inter['af_y'] >= 0.1) & (inter['af_y'] <= 0.9))]

df['color']=np.nan

# Load the example diamonds dataset

# Draw a scatter plot while assigning point colors and sizes to different
# variables in the dataset

colors = ["#ff82a0", "#d9c6e5", "#ffe0e6"]
sns.set_palette(sns.color_palette(colors))

custom_palette = {}

for q in set(df.dbsnp):
	afx = (np.average(df[df.dbsnp == q].af_x))
	afy = (np.average(df[df.dbsnp == q].af_y))

	if ((afx >= 0.1 and afx <= 0.9) and (afy <= 0.1 or afy >= 0.9)):
		df[df.dbsnp == q].color = 0
		custom_palette[q] = '#FF0B04'
	elif ((afy >= 0.1 and afy <= 0.9) and (afx <= 0.1 or afx >= 0.9)):
		df[df.dbsnp == q].color = 1
		custom_palette[q] = '#FF0B04'
	elif ((afy >= 0.1 and afy <= 0.9) and (afx >= 0.1 and afx <= 0.9)):
		df[df.dbsnp == q].color = 2
		custom_palette[q] = '#FF0B04'


for index in df.index:
	afx=df.loc[index, "af_x"]
	afy=df.loc[index, "af_y"]
	if ((afx >= 0.1 and afx <= 0.9) and (afy <= 0.1 or afy >= 0.9)):
		df.loc[index, "color"] = 0
	elif ((afy >= 0.1 and afy <= 0.9) and (afx <= 0.1 or afx >= 0.9)):
		df.loc[index, "color"] = 1
	elif ((afy >= 0.1 and afy <= 0.9) and (afx >= 0.1 and afx <= 0.9)):
		df.loc[index, "color"] = 2


sns.relplot(data= df, x="af_x", y="af_y", hue="color", legend=False)

plt.show()

#0.1 =< af_x =< 0.9 && (af_y =< 0.1 || af_y >= 0.9) = rosso

#0.1 =< af_y =< 0.9 && (af_x =< 0.1 || af_x >= 0.9) = blu

#0.1 =< af_x =< 0.9 && 0.1 =< af_y =< 0.9 = giallo

fil1 = "cov_filtering/notre/mock_TOT_1_S15.final.PILEUP.ASEQ.10dbsnp"
fil2 = "cov_filtering/notre/mock_TOT_3_S14.final.PILEUP.ASEQ.10dbsnp"
fil3 = "cov_filtering/notre/mock_TOT_6_S13.final.PILEUP.ASEQ.10dbsnp"

df1 = pd.read_csv(fil1, sep='\t')
df2 = pd.read_csv(fil2, sep='\t')
df3 = pd.read_csv(fil3, sep='\t')

set1 = set(df1.dbsnp.ravel())


set2 = set(df2.dbsnp.ravel())


set3 = set(df3.dbsnp.ravel())


sets = [set1, set2, set3]
labels = [fil1, fil2, fil3]
supervenn(sets, labels)

plt.show()

'''


colnames = ['val', 'dbsnp', 'boh']
df1 = pd.read_csv('p1', names=colnames, sep='\t')
df2 = pd.read_csv('p2', names=colnames, sep='\t')
df3 = pd.read_csv('p3', names=colnames, sep='\t')
df4 = pd.read_csv('p4', names=colnames, sep='\t')
#df1 = pd.read_csv('prova1', names=colnames, sep='\t')
#df2 = pd.read_csv('prova2', names=colnames, sep='\t')
#df3 = pd.read_csv('prova3', names=colnames, sep='\t')


dfs=[df1, df2, df3, df4]
df_x = df1.rename({'dbsnp' : 'rsid'}, axis=1)
df_y = df2.rename(columns = {'dbsnp' : 'rsid'})

result = pd.merge(df_x, df_y, on='rsid', how='outer', suffixes=('_x', '_y'))
print (result)
#result.to_csv("csv_prova", sep='\t')
'''
result1 = pd.merge(df4, df3, on='dbsnp', how='outer', suffixes=('_1', '_2'))
print (result1)
#result1.to_csv("csv_prova1", sep='\t')

result2 = pd.merge(result, result1, on='dbsnp', how='outer', suffixes=('_3', '_4'))
print (result2)
#result1.to_csv("csv_prova1", sep='\t')

'''

