import glob
from supervenn import supervenn
import matplotlib.pyplot as plt
from upsetplot import from_memberships,plot,from_contents
import itertools

out_pac = "/GROUPS/sharedRL/tirocinanti/elisa/out_pacbam"

def buildSets(globPath, minCov, minAF, maxAF, extention):
	sets = {}
	for f in glob.glob(globPath):
		fname = f.split("/")[2].split(".")[0].split("_")[0:3]
		fname = "_".join(fname)
		print(fname)
		sets[fname] = set()
		with open(f) as pileup:
			next(pileup)
			for l in pileup:
				l = l.strip().split()
				rsid = l[2]
				cov = int(l[10])
				af = float(l[9])
				if extention == "pabs" :
					rsid = l [1]
					cov = int(l[9])
					af = int(l[8])
	
				if cov >= minCov and af >= minAF and af <= maxAF:
					sets[fname].add(rsid)
	return sets

def plotUpset(data,outFile):
	data = from_contents(data)
	plot(data,sort_by="cardinality")
	plt.savefig(outFile)

def unionSubset(subset):
	sets = [subset[k] for k in subset]
	return set.intersection(*sets)

def unionSubsetCarbon(subset):
	sets = [subset[k] for k in subset if k != "Carbon_TOT_2"]
	return set.intersection(*sets)

minCov = [10,20]
af = [(0.1,0.9),(0.2,0.8)]
ex = [snps, pabs]

for cov, afRange, extention in itertools.product(minCov,af, ex):
	afMin = afRange[0]
	afMax = afRange[1]	
	setsCarbon = buildSets(out_pac+"/Carbon*."+extention,cov,afMin,afMax,extention)
	plotUpset(setsCarbon,f'./images/Carbon_{cov}_{afMin}_{afMax}_{extention}.png')
	setsMock = buildSets(out_pac+"/mock*."+extention,cov,afMin,afMax,extention)
	plotUpset(setsMock,f'./images/mock_{cov}_{afMin}_{afMax}_{extention}.png')
	setsXray = buildSets(out_pac+"/X-ray*."+extention,cov,afMin,afMax,extention)
	plotUpset(setsXray,f'./images/X-ray_{cov}_{afMin}_{afMax}_{extention}.png')
	setsProton = buildSets(out_pac+"/Proton*"+extention,cov,afMin,afMax,extention)
	plotUpset(setsProton,f'./images/Proton_{cov}_{afMin}_{afMax}_{extention}.png')	
	setTotal = {}
	setTotal["Carbon"] = unionSubset(setsCarbon)
	setTotal["Proton"] = unionSubset(setsProton)
	setTotal["X-ray"] = unionSubset(setsXray)
	setTotal["mock"] = unionSubset(setsMock)	
	plotUpset(setTotal, f'./images/allConditions_{cov}_{afMin}_{afMax}_{extention}.png')	
	setTotal["Carbon"] = unionSubsetCarbon(setsCarbon)
	plotUpset(setTotal, f'./images/allConditions_CarbonTOT2Removed_{cov}_{afMin}_{afMax}_{extention}.png')