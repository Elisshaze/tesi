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

fraction = ["POL", "TOT"]
in_pac = "/media/elisa/backup/out_pacbam/"
condition = ["mock", "Carbon", "Proton", "X-ray"]
minCov = 20
af = [(0.1,0.9),(0.2,0.8)]


def makeUnion (dfs) :
	merge = pd.merge(dfs[0], dfs[1], on='pos', how='outer', suffixes=('_0', '_1'))
	i=2
	for i in range(len(dfs)) :
		j= i+1
		result = pd.merge(dfs[i], dfs[i+1], on='pos', how='outer', suffixes=(f'_{str(i)}', f'_{str(j)}'))
		result1 = pd.merge(merge, result, n='pos', how='outer')
		merge = result1
	print (merge)




for con, frac, afRange in itertools.product(condition, fraction, af):
	afMin = afRange[0]
	afMax = afRange[1]
	path = f'{in_pac}{con}_{frac}*'
	print (path)

	dfs = []
	for f in glob.glob(path):
		print (f)
		df = pd.read_csv(f, sep='\t')
		first_column = df.pop('pos')
		df.insert(0, 'pos', first_column)
		dfs.append(df)
	
	makeUnion(dfs)