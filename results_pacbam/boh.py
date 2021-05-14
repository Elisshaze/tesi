import glob
from supervenn import supervenn
import matplotlib.pyplot as plt
from upsetplot import from_memberships,plot,from_contents
import itertools
import pandas as pd 
import os 
import seaborn as sns
import numpy as np


in_pac = "/media/elisa/backup/out_pacbam/"
out_scatter = "./snps_and_pabs/scatters/"
stats = out_scatter+"/stats"
out_venns = "./snps_and_pabs/venns/"


def buildSets(globPath, minCov, minAF, maxAF):
	sets = []
	names = []
	for f in glob.glob(globPath):
		extention = f.split(".")[2]
		fname = f.split("/")[5].split(".")[0].split("_")[0:3]
		fname = "_".join(fname)+"_"+extention
		print(fname)
		names.append(fname)
		with open(f) as pileup:
			next(pileup)
			subset = []
			for l in pileup:
				l = l.strip().split()
				rsid = l[1]
				if extention == "pabs" :
					cov = int(l[9])
					af = float(l[8])
				if extention == "snps" :
					cov = int(l[10])
					af = float(l[9])
				if cov >= minCov and af >= minAF and af <= maxAF:
					subset.append(rsid)
			s=set(subset)
			sets.append(s)
	return [sets, names]

colnames_snps = ['chr','pos','rsid','ref','alt','A','C','G','T','af','cov']
colnames_pabs= ['chr','pos','ref','alt','A','C','G','T','af','cov']

def makeScatter (snp, pab, name, cov, afMin, afMax):
	print ("in function make scatter")

	print (cov, afMin, afMax)
	df_snp = pd.read_csv(snp, sep='\t')
	df_snp = df_snp[(df_snp['cov'] >= cov)]
	print (df_snp)
	df_pab = pd.read_csv(pab, sep='\t')
	df_pab = df_pab[(df_pab['cov'] >= cov)]
	print (df_pab)

	inter = pd.merge(df_snp, df_pab, on='pos', how='inner')
	filt=((inter['af_x'] >= afMin) & (inter['af_x'] <= afMax)) | ((inter['af_y'] >= afMin) & (inter['af_y'] <= afMax))
	inter = inter[filt]
	print (inter)

	f = open (stats, 'a')
	lines = str(len(inter))
	f.write(name+'\t'+lines+'\n')
	f.close()

	#boh = sns.relplot(data= inter, x="af_x", y="af_y", s=5)
	boh=sns.jointplot(data= inter, x="af_x", y="af_y", s=7)
	#boh.set(xlabel="snps", ylabel="pabs")
	#plt.title(f'{name}_{afMin}_{afMax}_{cov}, {lines}' )
	plt.savefig(f'{out_scatter}{name}_{afMin}_{afMax}_{cov}.pdf', dpi=300)
	#plt.show()
	plt.close()


def makeVenn (sets, names, parameters):
	supervenn(sets, names, min_width_for_annotation=900, widths_minmax_ratio=0.1,
				col_annotations_area_height=1.2, rotate_col_annotations=True, sets_ordering='size')
	plt.title(parameters)
	#plt.show()
	plt.savefig(f'{out_venns}/{parameters}.pdf')
	plt.close()


condition = ["mock", "Carbon", "Proton", "X-ray"]
minCov = [50,80]
af = [(0.1,0.9),(0.2,0.8)]
ex = ['snps', 'pabs']
fraction = ['POL', 'TOT']
'''
for con, cov, afRange, extention in itertools.product(condition, minCov,af,ex):
	afMin = afRange[0]
	afMax = afRange[1]	
	
	for snp in glob.glob(f'{in_pac}{con}_*.snps'):
		pab = snp.replace('snps', 'pabs')
		
		fname = snp.split("/")[5].split(".")[0].split("_")[0:3]
		fname = "_".join(fname)
		#print (fname)
		#makeScatter(snp, pab, fname, cov, afMin, afMax)
'''	
for con, frac, cov, afRange,  in itertools.product(condition, fraction, minCov,af):
	afMin = afRange[0]
	afMax = afRange[1]

	pabs=[]
	print (f'{in_pac}{con}_{frac}*.pabs')
	for pab in glob.glob(f'{in_pac}{con}_{frac}_*.pabs'):
		
		pabs.append(pab)
	print (pabs)
	for j in range(3) :
		name1=pabs[j % len(pabs)].split("/")[5].split(".")[0].split("_")[0:3]
		name1="_".join(name1)
		
		name2=pabs[(j+1) % len(pabs)].split("/")[5].split(".")[0].split("_")[0:3]
		name2="_".join(name2)
		print(name1, name2)
		fname=name1+"*"+name2
		print (fname)
		makeScatter(pabs[j % len(pabs)], pabs[(j+1) % len(pabs)], fname, cov, afMin, afMax)

'''


for con, cov, afRange, in itertools.product(condition, minCov,af):
	afMin = afRange[0]
	afMax = afRange[1]
	parameters = f'{con}_{cov}_{afMin}_{afMax}'

	sets = buildSets(f'{in_pac}/{con}*',cov,afMin,afMax)
	setsCarbon = sets[0]
	namesCarbon = sets[1]

	makeVenn(setsCarbon, namesCarbon, parameters)
'''