import pandas as pd 
import glob
import os 
import plotly.graph_objects as go
from itertools import cycle
from supervenn import supervenn
import seaborn as sns
import matplotlib.pyplot as plt 
import numpy as np

pd.options.mode.chained_assignment = None  # default='warn'

def make_scatter(df,  condition, x, y, filt) :
	print ("in function make scatter")
	name = "results_cov/pdfs/"+condition+"_"+x+"_"+"_"+y+"_"+filt
	print (df, condition, x, y, filt)
	df['color']=np.nan

	colors = ["#ff82a0", "#d9c6e5", "#ffe0e6"]
	sns.set_palette(sns.color_palette(colors))

	for index in df.index:
		afx=df.loc[index, "af_x"]
		afy=df.loc[index, "af_y"]
		if ((afx >= 0.1 and afx <= 0.9) and (afy <= 0.1 or afy >= 0.9)):
			df.loc[index, "color"] = 0
		elif ((afy >= 0.1 and afy <= 0.9) and (afx <= 0.1 or afx >= 0.9)):
			df.loc[index, "color"] = 1
		elif ((afy >= 0.1 and afy <= 0.9) and (afx >= 0.1 and afx <= 0.9)):
			df.loc[index, "color"] = 2

	boh = sns.relplot(data= df, x="af_x", y="af_y", hue="color", legend=False)
	boh.set(xlabel=x, ylabel=y)
	plt.savefig(name+'.pdf', dpi=300)


def findnames (arr) :
	res=[]
	for fil in arr :
		name = fil.split("/")
		supp = name[2].split(".")
		res.append(supp[0])
	return res


mocks_n_cov10=[] #file mock notre, af non filtrato
mocks_n_cov20=[]
mocks_r_cov10=[]
mocks_r_cov20=[]

arrays=[mocks_n_cov10, mocks_r_cov10, mocks_n_cov20, mocks_r_cov20]
names_arrays=['mocks_n_cov10', 'mocks_r_cov10', 'mocks_n_cov20', 'mocks_r_cov20']
n_arrays=4


for file in glob.glob('cov_filtering/notre/mock_TOT*.10filtered'):
    mocks_n_cov10.append(file)

for file in glob.glob('cov_filtering/notre/mock_TOT*.20filtered'):
	mocks_n_cov20.append(file)

for file in glob.glob('cov_filtering/ricalib/mock_TOT*.10filtered'):
	mocks_r_cov10.append(file)

for file in glob.glob('cov_filtering/ricalib/mock_TOT*.20filtered'):
	mocks_r_cov20.append(file)

for array in arrays :
	array.sort()

for i in range(n_arrays) :
	curr=arrays[i] #determina un quale array, condizione, sto lavorando
	print curr
	names_curr = findnames(curr)
	print names_curr

	stats='results_cov/stats'+names_arrays[i]
	open (stats, 'w')
	cov='10' #quale coverage
	opt='n' #non ricalibrato (n), ricalibrato (r)
	if (i>1) :
		cov='20'
	if (i%2==1) :
		opt='r'

	dfs=[]
	for fil in  curr:
		df = pd.read_csv(fil, sep='\t')
		dfs.append(df)
	
	for j in range(3) :
		x=(j % len(dfs))
		y=((j+1) % len(dfs))

		df1=(dfs[j % len(dfs)])
		df2=(dfs[(j+1) % len(dfs)])
		
		inter = pd.merge(df1, df2, on='dbsnp', how='inner')

		#options for filtering on converage
		filA=((inter['af_x'] >= 0.1) & (inter['af_x'] <= 0.9)) | ((inter['af_y'] >= 0.1) & (inter['af_y'] <= 0.9))
		filB=((inter['af_x'] >= 0.2) & (inter['af_x'] <= 0.8)) | ((inter['af_y'] >= 0.2) & (inter['af_y'] <= 0.8))

		interA = inter[filA]
		make_scatter(interA, names_arrays[i], names_curr[x], names_curr[y], "A")
		interB = inter[filB]
		make_scatter(interA, names_arrays[i], names_curr[x], names_curr[y], "B")

	dfs=[]
		






