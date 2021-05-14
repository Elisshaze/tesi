import pandas as pd 
from pandas import *
import glob
import os 
import plotly.graph_objects as go
from supervenn import supervenn
import matplotlib.pyplot as plt 
import re
import numpy as np
import seaborn as sns
import itertools

out_pacbam = "/media/elisa/backup/out_pacbam/"
out_aseq = "/media/elisa/backup/out_aseq/"
out_aseq_recalib = "/media/elisa/backup/out_aseq_recalib/"

minCov = 20
def makeDfs (files, names, indexes, afMin, afMax, top):
	dfs = []
	dfs_filt = []
	intersect_names = []
	

	for f in files :
		df = pd.read_csv(f, sep='\t')
		df = df[(df['cov'] >= minCov)]
		dfs.append(df)
		#filt=((df['af'] >= afMin) & (df['af'] <= afMax)) & (df['cov'] >= minCov)
		#df = df[filt]
		#dfs_filt.append(df)
	res = []
	

	s = top #second
	intersections = []
	names_couples = []
	for i in range(1, 4):
		print (i)
		f = s # first
		s = (top - (i%3))
		fi = "S"+str(f)
		se = "S"+str(s)
		print (fi, se)
		#files[indexes.index()]

		#aseq * aseq
		df = pd.merge(dfs[indexes.index(fi+".ASEQ")], dfs[indexes.index(se+".ASEQ")], on="dbsnp", how="inner")
		filt=((df['af_x'] >= afMin) & (df['af_x'] <= afMax)) & ((df['af_y'] >= afMin) & (df['af_y'] <= afMax))
		df = df[filt]
		print(df)
		intersections.append(df)
		names_couples.append([names[indexes.index(fi+".ASEQ")], names[indexes.index(se+".ASEQ")]])

		#aseq * aseq (almeno uno dei due)
		inter = pd.merge(dfs[indexes.index(fi+".ASEQ")], dfs[indexes.index(se+".ASEQ")], on='dbsnp', how='inner')
		filt=((inter['af_x'] >= afMin) & (inter['af_x'] <= afMax)) | ((inter['af_y'] >= afMin) & (inter['af_y'] <= afMax))
		inter = inter[filt]
		print(inter)
		intersections.append(df)
		names_couples.append([names[indexes.index(fi+".ASEQ")], names[indexes.index(se+".ASEQ")]])
		#recalib * aseq
		inter = pd.merge(dfs[indexes.index(fi+".recalibrated")], dfs[indexes.index(se+".ASEQ")], on="dbsnp", how="inner")
		inter = inter[filt]
		print(inter)
		intersections.append(df)
		names_couples.append([names[indexes.index(fi+".recalibrated")], names[indexes.index(se+".ASEQ")]])

		#aseq * recalib
		inter= pd.merge(dfs[indexes.index(fi+".ASEQ")], dfs[indexes.index(se+".recalibrated")], on="dbsnp", how="inner")
		inter = inter[filt]
		print(inter)
		intersections.append(df)
		names_couples.append([names[indexes.index(fi+".ASEQ")], names[indexes.index(se+".recalibrated")]])

		#snps * snps
		inter = pd.merge(dfs[indexes.index(fi+".snps")], dfs[indexes.index(se+".snps")],  on="rsid", how="inner")
		inter = inter[filt]
		print(inter)
		intersections.append(df)
		names_couples.append([names[indexes.index(fi+".snps")], names[indexes.index(se+".snps")]])

		#pabs * pabs
		inter = pd.merge(dfs[indexes.index(fi+".pabs")], dfs[indexes.index(se+".pabs")], on="pos", how="inner")
		inter = inter[filt]
		print(inter)
		intersections.append(df)
		names_couples.append([names[indexes.index(fi+".pabs")], names[indexes.index(se+".pabs")]])

		#1_aseq * 1_snps
		df_renamed = dfs[indexes.index(fi+".ASEQ")].rename(columns = {'dbsnp' : 'rsid'}, inplace = False)
		inter = pd.merge(df_renamed, dfs[indexes.index(fi+".snps")], on="rsid", how="inner")
		inter = inter[filt]
		print(inter)
		intersections.append(df)
		names_couples.append([names[indexes.index(fi+".ASEQ")], names[indexes.index(fi+".snps")]])

		#2_aseq * 2_snps
		df_renamed= dfs[indexes.index(se+".ASEQ")].rename(columns = {'dbsnp' : 'rsid'}, inplace = False)
		inter = pd.merge(df_renamed, dfs[indexes.index(se+".snps")], on="rsid", how="inner")
		inter = inter[filt]
		print(inter)
		intersections.append(df)
		names_couples.append([names[indexes.index(se+".ASEQ")], names[indexes.index(se+".snps")]])

		print (DataFrame(names_couples))
		intersections = []
		
		names_couples = []
		


def makeScatters (names, dfss, afMin, afMax)	:
	for dfs in dfss:
		df1 = dfs[0]
		df2 = dfs[1]
		inter = pd.merge(df1, df2, on='pos', how='inner')
		filt=((inter['af_x'] >= afMin) & (inter['af_x'] <= afMax)) | ((inter['af_y'] >= afMin) & (inter['af_y'] <= afMax))
		inter = inter[filt]
		print (inter)

		
'''
mock_tot 1 mock_tot3
mocktot1.aseq * mocktot3.aseq 		mocktot1.recalib * mocktot3.recalib    mock_tot1.recalib* mock_tot3.aseq 
mocktot1.aseq* mocktot3.recalib		1*3, almeno uno compreso				mock_tot1.pabs * mock_tot3.pabs
mocktot1.snps * mocktot3.snps		mocktot1.aseq * mocktot1.snps			mock_tot3.snps * mock_tot3.aseq
'''

def retrieveName (files):
	names =[]
	indexes = []
	for f in files:
		split = f.split(".")
		extention = split[-1]
		control = split[-3]
		fname = f.split("/")[5].split(".")[0].split("_")[0:3]
		index = f.split("/")[5].split(".")[0].split("_")[3]+"."+extention
		
		base= "_".join(fname)
		fname = base+"_"+extention
		if control == "recalibrated" :
			fname = base+"_"+control+"."+extention
			index = f.split("/")[5].split(".")[0].split("_")[3]+"."+control
		names.append(fname)
		indexes.append(index)	
	return names, indexes


paths = [out_aseq, out_aseq_recalib, out_pacbam]
#condition = ["mock", "Carbon", "Proton", "X-ray"]
condition = ["mock", "Carbon"]
#fraction = ["POL", "TOT"]
fraction = ["TOT"]
#af = [(0.1,0.9),(0.2,0.8)]
af = [(0.1,0.9)]

for con, frac, afRange in itertools.product(condition, fraction, af):
	afMin = afRange[0]
	afMax = afRange[1]
	fil = f'{con}_{frac}'

	files = []
	for path in paths :
		path = f'{path}{fil}_*'
		for f in glob.glob(path):
			#print (f)
			files.append(f)
			'''
			df = pd.read_csv(f, sep='\t')
			filt=((df['af'] >= afMin) & (df['af'] <= afMax)) & (df['cov'] >= minCov)
			df = df[filt]
			first_column = df.pop('pos')
			df.insert(0, 'pos', first_column)
			dfs.append(df)
			'''
		#makeUnion(dfs, f'{files}_{afMin}_{afMax}')
	
	files.sort()
	names, indexes = retrieveName(files)
	print (indexes)
	print (int(indexes[0].split(".")[0].split("S")[1]))
	#print(files)
	makeDfs(files, names, indexes, afMin, afMax, (int(indexes[0].split(".")[0].split("S")[1])))