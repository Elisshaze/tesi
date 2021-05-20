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
pd.set_option('display.max_colwidth', None)
out_pacbam = "/media/elisa/backup/out_pacbam/"
out_aseq = "/media/elisa/backup/out_aseq/"
out_aseq_recalib = "/media/elisa/backup/out_aseq_recalib/"



#((afy >= afMin and afy <= afMax) and (afx >= afMin and afx <= afMax))


minCov = 20
def makeDfs (files, names, indexes, con, afMin, afMax, top):
	dfs = []
	prova = []
	i=0
	for f in files :
		df = pd.read_csv(f, sep='\t')
		df = df[(df['cov'] >= minCov)]
		dfs.append(df)
		prova.append(f)
		#filt=((df['af'] >= afMin) & (df['af'] <= afMax)) & (df['cov'] >= minCov)
		#df = df[filt]
		#dfs_filt.append(df)
		i+=1

	print (DataFrame(prova))
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

	
		#aseq * aseq
		df = pd.merge(dfs[indexes.index(fi+".ASEQ")], dfs[indexes.index(se+".ASEQ")], on="dbsnp", how="inner")
		filt=((df['af_x'] >= afMin) & (df['af_x'] <= afMax)) & ((df['af_y'] >= afMin) & (df['af_y'] <= afMax))
		df = df[filt]
		intersections.append(df)
		names_couples.append([names[indexes.index(fi+".ASEQ")], names[indexes.index(se+".ASEQ")]])

		#aseq * aseq (almeno uno dei due)
		inter = pd.merge(dfs[indexes.index(fi+".ASEQ")], dfs[indexes.index(se+".ASEQ")], on='dbsnp', how='inner')
		filt=((inter['af_x'] >= afMin) & (inter['af_x'] <= afMax)) | ((inter['af_y'] >= afMin) & (inter['af_y'] <= afMax))
		inter = inter[filt]
		intersections.append(inter)
		names_couples.append([names[indexes.index(fi+".ASEQ")], names[indexes.index(se+".ASEQ")]])

		#recalib * aseq
		inter = pd.merge(dfs[indexes.index(fi+".recalibrated")], dfs[indexes.index(se+".ASEQ")], on="dbsnp", how="inner")
		filt=((inter['af_x'] >= afMin) & (inter['af_x'] <= afMax)) | ((inter['af_y'] >= afMin) & (inter['af_y'] <= afMax))
		inter = inter[filt]
		intersections.append(inter)
		names_couples.append([names[indexes.index(fi+".recalibrated")], names[indexes.index(se+".ASEQ")]])

		#aseq * recalib
		inter= pd.merge(dfs[indexes.index(fi+".ASEQ")], dfs[indexes.index(se+".recalibrated")], on="dbsnp", how="inner")
		filt=((inter['af_x'] >= afMin) & (inter['af_x'] <= afMax)) | ((inter['af_y'] >= afMin) & (inter['af_y'] <= afMax))
		inter = inter[filt]
		intersections.append(inter)
		names_couples.append([names[indexes.index(fi+".ASEQ")], names[indexes.index(se+".recalibrated")]])

		#snps * snps
		inter = pd.merge(dfs[indexes.index(fi+".snps")], dfs[indexes.index(se+".snps")],  on="rsid", how="inner")
		filt=((inter['af_x'] >= afMin) & (inter['af_x'] <= afMax)) | ((inter['af_y'] >= afMin) & (inter['af_y'] <= afMax))
		inter = inter[filt]
		intersections.append(inter)
		names_couples.append([names[indexes.index(fi+".snps")], names[indexes.index(se+".snps")]])

		#pabs * pabs
		inter = pd.merge(dfs[indexes.index(fi+".pabs")], dfs[indexes.index(se+".pabs")], on="pos", how="inner")
		filt=((inter['af_x'] >= afMin) & (inter['af_x'] <= afMax)) | ((inter['af_y'] >= afMin) & (inter['af_y'] <= afMax))
		inter = inter[filt]
		intersections.append(inter)
		names_couples.append([names[indexes.index(fi+".pabs")], names[indexes.index(se+".pabs")]])
		
		#1_aseq * 1_snps
		df_renamed = dfs[indexes.index(fi+".ASEQ")].rename(columns = {'dbsnp' : 'rsid'})
		inter = pd.merge(df_renamed, dfs[indexes.index(fi+".snps")], on="rsid", how="inner")		
		filt=((inter['af_x'] >= afMin) & (inter['af_x'] <= afMax)) | ((inter['af_y'] >= afMin) & (inter['af_y'] <= afMax))
		inter = inter[filt]
		intersections.append(inter)
		names_couples.append([names[indexes.index(fi+".ASEQ")], names[indexes.index(fi+".snps")]])

		#2_aseq * 2_snps
		df_renamed= dfs[indexes.index(se+".ASEQ")].rename(columns = {'dbsnp' : 'rsid'})
		inter = pd.merge(df_renamed, dfs[indexes.index(se+".snps")], on="rsid", how="inner")		
		filt=((inter['af_x'] >= afMin) & (inter['af_x'] <= afMax)) | ((inter['af_y'] >= afMin) & (inter['af_y'] <= afMax))
		inter = inter[filt]
		intersections.append(inter)
		names_couples.append([names[indexes.index(se+".ASEQ")], names[indexes.index(se+".snps")]])

		#print(intersections)
		title = f'{con}_{afMin}-{afMax}'
		makeScatters(intersections, names_couples, afMin, afMax, title)
		MakeKDE(intersections, names_couples, afMax, title)
		intersections = []
		names_couples = []
		

def makeScatters (dfs, names_couples, afMin, afMax, title) :
	fig = plt.figure()
	fig.set_size_inches(30, 30)
	fig.subplots_adjust(hspace=0.3, wspace=0.3)
	i = 1
	for df in dfs:
		lines = str(len(df))
		ax = fig.add_subplot(4, 2, i)
		ax.autolayout=False
		ax.set_ylim([-0.1, 1.2])
		if i== 1 : 
			ax.set_title(names_couples[i-1][0]+"*"+names_couples[i-1][1]+" both afs in range, "+lines, fontsize=24, ha='center')
		else :
			ax.set_title(names_couples[i-1][0]+"*"+names_couples[i-1][1]+", "+lines, fontsize=24, ha='center')
		
	
		df['color']=np.nan

		colors = ["#13f007", "#ebe156", "#4d32ed"]
		sns.set_palette(sns.color_palette(colors))
		for index in df.index:
			afx=df.loc[index, "af_x"]
			afy=df.loc[index, "af_y"]
			if ((afx >= afMin and afx <= afMax) and (afy < afMin or afy > afMax)):
				df.loc[index, "color"] = 0
			elif ((afy >= afMin and afy <= afMax) and (afx < afMin or afx > afMax)):
				df.loc[index, "color"] = 1
			else :
				df.loc[index, "color"] = 2
				
		sns.scatterplot(data = df,  x="af_x", y="af_y", hue="color",  s=40, legend=False, ax=ax)
		i+=1 
	plt.savefig(f"./results_pacbam/grid/{title}_scatter.pdf", dpi=1000)
	#plt.show()

def MakeKDE (dfs, names_couples, afMax, title) :
	fig = plt.figure()
	fig.set_size_inches(30, 30)
	fig.subplots_adjust(hspace=0.3, wspace=0.3)
	i = 1
	for df in dfs:
		lines = str(len(df))
		ax = fig.add_subplot(4, 2, i)
		ax.autolayout=False
		ax.set_ylim([-0.1, 2.5])
		if afMax == 0.9 :
			ax.set_ylim([-0.1, 3.7])
		if i== 1 : 
			ax.set_title(names_couples[i-1][0]+"*"+names_couples[i-1][1]+" both afs in range, "+lines, fontsize=24, ha='center')
		else :
			ax.set_title(names_couples[i-1][0]+"*"+names_couples[i-1][1]+", "+lines, fontsize=24, ha='center')
		df1=df[['af_x', 'af_y']]
		sns.kdeplot(data = df1, bw_adjust=.3, ax=ax)
		i+=1 
	plt.savefig(f"./results_pacbam/grid/{title}_kde.pdf", dpi=1000)
	#plt.show()
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
condition = ["mock", "Carbon", "Proton", "X-ray"]
#condition = ["mock", "Carbon"]
fraction = ["POL", "TOT"]
af = [(0.1,0.9),(0.2,0.8)]
#af = [(0.1,0.9)]

for con, frac, afRange in itertools.product(condition, fraction, af):
	afMin = afRange[0]
	afMax = afRange[1]
	fil = f'{con}_{frac}'

	files = []
	for path in paths :
		path = f'{path}{fil}_*'
		for f in glob.glob(path):
			files.append(f)
	
	files.sort()
	print (files, )
	names, indexes = retrieveName(files)
	print (names, indexes)
	print (afMin, afMax)
	makeDfs(files, names, indexes, fil, afMin, afMax, (int(indexes[0].split(".")[0].split("S")[1])))