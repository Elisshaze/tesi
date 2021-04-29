import pandas as pd 
import glob
import os 
import plotly.graph_objects as go
from itertools import cycle
from supervenn import supervenn
import matplotlib.pyplot as plt 

##############################################################################################
#creazione dei diagrammi di venn usando supervenn
#utilizzo i file contenenti solo il campo -dbsnp-, per cov >= 10 e cov >= 20
#No dataframe ma sets


##############################################################################################

colnames = ['chr','pos','dbsnp','MAF','ref','alt','A','C','G','T','RD','Ars','Crs','Grs','Trs','af','cov']

mocks_n_cov10=[] #file mock notre, af non filtrato
mocks_n_cov20=[]
mocks_r_cov10=[]
mocks_r_cov20=[]

arrays=[mocks_n_cov10, mocks_r_cov10, mocks_n_cov20, mocks_r_cov20]
names_arrays=['mocks_n_cov10', 'mocks_r_cov10', 'mocks_n_cov20', 'mocks_r_cov20']
n_arrays=4


#options for filtering on converage
#filA=((intersect['af_x'] >= 0.1) & (intersect['af_x'] <= 0.9)) | ((intersect['af_y'] >= 0.1) & (intersect['af_y'] <= 0.9))
#filB=((intersect['af_x'] >= 0.2) & (intersect['af_x'] <= 0.8)) | ((intersect['af_y'] >= 0.2) & (intersect['af_y'] <= 0.8))


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

print arrays


for i in range(n_arrays) :
	curr=arrays[i] #determina un quale array, condizione, sto lavorando
	
	stats='results_cov/stats'+names_arrays[i]
	open (stats, 'w')
	cov='10' #quale coverage
	opt='n' #non ricalibrato (n), ricalibrato (r)
	if (i>1) :
		cov='20'
	if (i%2==1) :
		opt='r'

	df = 'df_'+ opt + cov
	print df
	dfs=[]
	for fil in arrays[1] :
		df = pd.read_csv(fil, names=colnames, sep='\t', low_memory=False )
		print df







