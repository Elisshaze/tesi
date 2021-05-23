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

out_pacbam = "/media/elisa/backup/out_pacbam_recalib/"
out_aseq = "/media/elisa/backup/out_aseq/"
out_aseq_recalib = "/media/elisa/backup/out_aseq_recalib/"
out_grid_wes = "./grid_WES2"
out_wes = "./out_pacbam_WES"

minCov = 20 
filesWES = [] #array containing path to SRR files
dfsWES = [] #array containg the SRR filtered

def processSRR(afMin, afMax):
	
	SRR98pabs = "./out_pacbam_WES/SRR9166098.sorted.groups.dedup.realigned.recalibrated.pabs"
	SRR98snps = "./out_pacbam_WES/SRR9166098.sorted.groups.dedup.realigned.recalibrated.snps"

	SRR99pabs = "./out_pacbam_WES/SRR9166099.sorted.groups.dedup.realigned.recalibrated.pabs"
	SRR99snps = "./out_pacbam_WES/SRR9166099.sorted.groups.dedup.realigned.recalibrated.snps"
	filesWES = [SRR98pabs, SRR98snps, SRR99pabs, SRR99snps]

	df98p = pd.read_csv(SRR98pabs, sep='\t')
	df98s = pd.read_csv(SRR98snps, sep='\t')
	df99p = pd.read_csv(SRR99pabs, sep='\t')
	df99s = pd.read_csv(SRR99snps, sep='\t')

	dfsWES=[df98p, df98s, df99p, df99s]
	return filesWES, dfsWES
'''
	inter = pd.merge(df98p, df99p, on= "pos", how="inner")
	filt=((inter['af_x'] >= afMin) & (inter['af_x'] <= afMax)) | ((inter['af_y'] >= afMin) & (inter['af_y'] <= afMax)) & ((inter['cov_x'] >= minCov) & (inter['cov_y'] >= minCov))
	df98p99p = inter[filt]
	dfsWES.append(df98p99p)
	print (df98p99p)
'''
	



####################################################################################
def retrieveName (files, index):
	names =[]
	for f in files:
		split = f.split(".")
		extention = split[-1]
		control = split[-3]
		fname = f.split("/")[index].split(".")[0].split("_")[0:3]

		base= "_".join(fname)
		fname = base+"."+extention
		if control == "recalibrated" :
			fname = base+"_"+control+"."+extention
		names.append(fname)
	return names

def MakeScattersWES(dfs, names, title):
	fig = plt.figure()
	fig.set_size_inches(30, 30)
	fig.subplots_adjust(hspace=0.3, wspace=0.3)
	i = 0
	j = 1
	for df in dfs:
		lines = str(len(df))
		ax = fig.add_subplot(2, 2, j)
		#ax.autolayout=False
		ax.set_ylim([-0.1, 1.2])
		print ("OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO", i)
		ax.set_title(names[i][0]+"*"+names[i][1]+", "+lines, fontsize=20, ha='center')
		df['color']=np.nan
		colors = ['#747FE3', '#8EE35D', '#E37346']
		for index in df.index:
			afx=df.loc[index, "af_x"]
			afy=df.loc[index, "af_y"]
			if ((afx >= afMin and afx <= afMax) and (afy < afMin or afy > afMax)):
				df.loc[index, "color"] = "only_X"
			elif ((afy >= afMin and afy <= afMax) and (afx < afMin or afx > afMax)):
				df.loc[index, "color"] = "only_Y"
			elif ((afy >= afMin and afy <= afMax) and (afx >= afMin or afx <= afMax)):
				df.loc[index, "color"] = "both"
				
		sns.scatterplot(data = df,  x="af_x", y="af_y", hue="color", palette = colors, s=17, legend=True, ax=ax)

		ax1 = fig.add_subplot(2, 2, j+1)
		#ax.autolayout=False
		ax1.set_ylim([-0.1, 2])
		ax1.set_title(names[i][0]+"*"+names[i][1]+", "+lines, fontsize=20, ha='center')
		df1=df[['af_x', 'af_y']]
		sns.kdeplot(data = df1, bw_adjust=.3, ax=ax1)


		i+=1
		j+=2 
	#plt.show()

	plt.savefig(f"{out_grid_wes}/{title}_scatter+kde.pdf", dpi=1000)

def StatsWes() : 
	cov=[10, 20]
	af = [(0.1,0.9),(0.2,0.8)]
	

	for afRange, cv in itertools.product(af, cov) :
		print (afRange, cv)
		minCov = cv
		afMin = afRange[0]
		afMax = afRange[1]

		couples = []
		dfs = []

		for i in range(2):
			inter = pd.merge(dfsWES[i], dfsWES[i+2], on="pos", how='inner')
			filt=((inter['af_x'] >= afMin) & (inter['af_x'] <= afMax)) | ((inter['af_y'] >= afMin) & (inter['af_y'] <= afMax)) & ((inter['cov_x'] >= minCov) & (inter['cov_y'] >= minCov))
			inter = inter[filt]
			dfs.append(inter)
			couple = [namesWES[i], namesWES[i+2]]
			couples.append(couple)
			
		title = f'WES_{minCov}_{afMin}-{afMax}'
		MakeScattersWES(dfs, couples, title)

def makePlots(index, mode, files, names, afMin, afMax, fil, dfsWES):
	print (dfsWES)
	pivot = dfsWES[index]
	print (pivot)
	print("INDEXXXXXXXXXXXXXXX", index, namesWES[index])
	name_pivot = namesWES[index]
	couples = []
	dfs = []

	for f in files:
		df = pd.read_csv(f, sep='\t')
		inter = pd.merge(pivot, df, on='pos', how='inner')
		filt=((inter['af_x'] >= afMin) & (inter['af_x'] <= afMax)) | ((inter['af_y'] >= afMin) & (inter['af_y'] <= afMax)) & ((inter['cov_x'] >= minCov) & (inter['cov_y'] >= minCov))
		inter = inter[filt]
		dfs.append(inter)
		couples.append([name_pivot, names[files.index(f)]])
	title = f'{name_pivot}*{fil}'
	makeScatters (dfs, couples, afMin, afMax, title)
	MakeKDE  (dfs, couples, afMax, title)
	
def makeScatters (dfs, names_couples, afMin, afMax, title) :
	fig = plt.figure()
	fig.set_size_inches(30, 30)
	fig.subplots_adjust(hspace=0.3, wspace=0.3)
	i = len(names_couples)
	for df in reversed(dfs):
		lines = str(len(df))
		ax = fig.add_subplot(5, 2, i)
		#ax = fig.add_subplot(1, 1, 1)
		ax.autolayout=False
		ax.set_ylim([-0.1, 1.2])
		ax.set_title(names_couples[i-1][0]+"*"+names_couples[i-1][1]+", "+lines, fontsize=30, ha='center')
		
		df['color']=np.nan
		
		for index in df.index:
			afx=df.loc[index, "af_x"]
			afy=df.loc[index, "af_y"]
			if ((afx >= afMin and afx <= afMax) and (afy < afMin or afy > afMax)):
				df.loc[index, "color"] = "only_X"
			elif ((afy >= afMin and afy <= afMax) and (afx < afMin or afx > afMax)):
				df.loc[index, "color"] = "only_Y"
			elif ((afy >= afMin and afy <= afMax) and (afx >= afMin or afx <= afMax)):
				df.loc[index, "color"] = "both"
				
		sns.scatterplot(data = df,  x="af_x", y="af_y", hue="color",  s=40, ax=ax)
		i-=1 
	plt.savefig(f"{out_grid_wes}/{title}_{afMax}_scatter.pdf", dpi=1000)
	plt.show()

def MakeKDE (dfs, names_couples, afMax, title) :
	fig = plt.figure()
	fig.set_size_inches(30, 30)
	fig.subplots_adjust(hspace=0.3, wspace=0.3)
	i = len(names_couples)
	for df in reversed(dfs):
		lines = str(len(df))
		ax = fig.add_subplot(5, 2, i)
		ax.autolayout=False
		#ax.set_ylim([-0.1, 6])
		#if afMax == 0.9 :
			#ax.set_ylim([-0.1, 6])
		
		ax.set_title(names_couples[i-1][0]+"*"+names_couples[i-1][1]+", "+lines, fontsize=24, ha='center')
		df1=df[['af_x', 'af_y']]
		sns.histplot(data = df1, bins = 30, kde = True, ax=ax)
		i-=1 
	plt.show()
	plt.savefig(f"{out_grid_wes}/{title}_{afMax}_kde.pdf", dpi=1000)

def buildSets(files, minCov, minAF, maxAF, mode):
	sets = []
	names = []
	for f in files:
		extention = f.split(".")[-1]
		fname = f.split("/")[5].split(".")[0].split("_")[0:3]
		fname = "_".join(fname)+"_"+extention
		print("BUILD SETSSSS", fname)
		names.append(fname)
		with open(f) as pileup:
			next(pileup)
			subset = []
			for l in pileup:
				l = l.strip().split()
				if mode == "p" :
					rsid = l[1]
					if extention == "pabs" :
						cov = int(l[9])
						af = float(l[8])
					if extention == "ASEQ" :
						cov = int(l[16])
						af = float(l[15])
				else :
					rsid=l[2]
					if extention == "snps" :
						cov = int(l[10])
						af = float(l[9])
					if extention == "ASEQ" :
						cov = int(l[16])
						af = float(l[15])
				if cov >= minCov and af >= minAF and af <= maxAF:
					subset.append(rsid)
			s=set(subset)
			sets.append(s)
	return [sets, names]

def makeVenn (sets, names, parameters):
	supervenn(sets, names, min_width_for_annotation=900, widths_minmax_ratio=0.1,
				col_annotations_area_height=1.2, rotate_col_annotations=True, sets_ordering='size')
	plt.title(parameters)
	plt.show()
	plt.savefig(f'/venns_wes/{parameters}.pdf')
	plt.close()


def Main () :
	condition = ["mock", "Carbon", "Proton", "X-ray"]
	#condition = ["mock"]
	fraction = ["POL", "TOT"]
	af = [(0.1,0.9),(0.2,0.8)]
	pac = ['p', 's', 'p', 's'] 
	#af = [(0.1,0.9)]

	i=0
	for  con, frac, afRange, p in itertools.product(condition, fraction, af, pac):
		files = []
		afMin = afRange[0]
		afMax = afRange[1]
		fil = f'{con}_{frac}'

		if p == "s" :
			for f in glob.glob(f'{out_pacbam}{fil}*.snps') :
				files.append(f)
		else :
			for f in glob.glob(f'{out_pacbam}{fil}*.pabs') :
				files.append(f)	

		for f in glob.glob(f'{out_aseq}{fil}*') :
				files.append(f)	
		for f in glob.glob(f'{out_aseq_recalib}{fil}*') :
				files.append(f)
		
		files.sort()
		print (files)
		names = retrieveName(files, 5)
		print (names)
		print (afMin, afMax)
		#makePlots
		#(files, names, indexes, fil, afMin, afMax, (int(indexes[0].split(".")[0].split("S")[1])))
		if afMin == 0.1 :
			dfsWES= dfs_WES1
		else :
			dfsWES = dfs_WES2
		makePlots(i%4, p, files, names, afMin, afMax, fil, dfsWES)
		#sets, names = buildSets(files, minCov, afMin, afMax, p)
		i+=1

#Retrieve and make dfs of the WES files
SRR98pabs = "./out_pacbam_WES/SRR9166098.sorted.groups.dedup.realigned.recalibrated.pabs"
SRR98snps = "./out_pacbam_WES/SRR9166098.sorted.groups.dedup.realigned.recalibrated.snps"

SRR99pabs = "./out_pacbam_WES/SRR9166099.sorted.groups.dedup.realigned.recalibrated.pabs"
SRR99snps = "./out_pacbam_WES/SRR9166099.sorted.groups.dedup.realigned.recalibrated.snps"
filesWES = [SRR98pabs, SRR98snps, SRR99pabs, SRR99snps]

df98p = pd.read_csv(SRR98pabs, sep='\t')
df98s = pd.read_csv(SRR98snps, sep='\t')
df99p = pd.read_csv(SRR99pabs, sep='\t')
df99s = pd.read_csv(SRR99snps, sep='\t')

dfsWES=[df98p, df98s, df99p, df99s]
namesWES = retrieveName(filesWES, 2)

af = [(0.1,0.9),(0.2,0.8)]
i=0
dfs_WES1=[]
dfs_WES2=[]

#for each df, filter
for afRange in af :
	afMin = afRange[0]
	afMax = afRange[1]
	print (afMin, afMax)

	for df in dfsWES :
		df = df[((df["af"]>= afMin) & (df["af"] <= afMax) & (df["cov"] >= 20))]

		if i < 4:
			dfs_WES1.append(df)
		else :
			dfs_WES2.append(df)
	i +=1


#DNA consensus: maintain only entries found in both DNA with af in range.
i=1
for afRange in (af):
	if i==1:
		a=dfs_WES1 # dfs with range 0.1-0.9
	else :
		a=dfs_WES2 #dfs with range 0.2-0.8
	afMin = afRange[0]
	afMax = afRange[1]
	print (afMin, afMax)
	#PABS
	inter = pd.merge(a[0], a[2], on= "pos", how="inner")
	print (inter) #intersection
	df98p_filt = df98p.pos.isin(inter.pos) #filtering for 98p
	df99p_filt = df99p.pos.isin(inter.pos) #fitlering for 99p
	df98p_new = df98p[df98p_filt]
	df99p_new = df99p[df99p_filt]

	print (df98p_filt, df99p_filt, df98p_new, df99p_new)

	inter = pd.merge(df98p, df98p_new, how="inner", on="pos")
	sns.scatterplot(data=inter, x="af_x", y="af_y")
	#questo scatterplot non mi convinve nemmeno un po'
	plt.show()

	#SNPS
	inter = pd.merge(a[1], a[3], on= "pos", how="inner")
	df98s_filt = df98s.pos.isin(inter.pos)
	df99s_filt = df99s.pos.isin(inter.pos)
	df98s_new = df98p[df98s_filt]
	df99s_new = df99p[df99s_filt]

	if i == 0:
		dfs_WES1=[df98p_new, df98s_new, df99p_new, df99s_new]
	else :
		dfs_WES2=[df98p_new, df98s_new, df99p_new, df99s_new]
	i+=1

Main()


