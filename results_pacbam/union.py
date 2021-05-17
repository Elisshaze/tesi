import pandas as pd 
import glob
import os 
import plotly.graph_objects as go
from supervenn import supervenn
import matplotlib.pyplot as plt 
import re
import numpy as np
import seaborn as sns
import itertools


in_pac = "/media/elisa/backup/out_pacbam/"
out_csvs = "./csvs/"
'''
def makeScatters ():
	fig = plt.figure()
	fig.subplots_adjust(hspace=0.4, wspace=0.4)
	for i in range(1, 7):
    	ax = fig.add_subplot(2, 3, i)
    	ax.text(0.5, 0.5, str((2, 3, i)),
           		fontsize=18, ha='center')
'''

def makeUnion (dfs, path) :
	merge = pd.merge(dfs[0], dfs[1], on='pos', how='outer', suffixes=('_0', '_1'))
	i=2
	while i < (len(dfs)-1) :
		j= i+1
		result = pd.merge(dfs[i], dfs[i+1], on='pos', how='outer', suffixes=(f'_{str(i)}', f'_{str(j)}'))
		result1 = pd.merge(merge, result, on='pos', how='outer')
		merge = result1
		i+=2
	print (merge)
	merge.to_csv(f'{out_csvs}{path}.csv')


condition = ["mock", "Carbon", "Proton", "X-ray"]
minCov = 20
af = [(0.1,0.9),(0.2,0.8)]
fraction = ["POL", "TOT"]

for con, frac, afRange in itertools.product(condition, fraction, af):
	afMin = afRange[0]
	afMax = afRange[1]
	files = f'{con}_{frac}'
	path = f'{in_pac}{files}_*'

	print (path)

	dfs = []
	for f in glob.glob(path):
		print (f)
		
		df = pd.read_csv(f, sep='\t')
		filt=((df['af'] >= afMin) & (df['af'] <= afMax)) & (df['cov'] >= minCov)
		df = df[filt]
		first_column = df.pop('pos')
		df.insert(0, 'pos', first_column)
		dfs.append(df)
		
	makeUnion(dfs, f'{files}_{afMin}_{afMax}')
