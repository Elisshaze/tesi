import glob
from supervenn import supervenn
import matplotlib.pyplot as plt
from upsetplot import from_memberships,plot,from_contents
import itertools

def buildSets(globPath, minCov, minAF, maxAF):
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
                cov = int(l[16])
                af = float(l[15])
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




for cov, afRange in itertools.product(minCov,af):
    afMin = afRange[0]
    afMax = afRange[1]

    setsCarbon = buildSets("../out_aseq/Carbon*.ASEQ",cov,afMin,afMax)
    plotUpset(setsCarbon,f'./images/Carbon_{cov}_{afMin}_{afMax}.png')
    setsMock = buildSets("../out_aseq/mock*.ASEQ",cov,afMin,afMax)
    plotUpset(setsMock,f'./images/mock_{cov}_{afMin}_{afMax}.png')
    setsXray = buildSets("../out_aseq/X-ray*.ASEQ",cov,afMin,afMax)
    plotUpset(setsXray,f'./images/X-ray_{cov}_{afMin}_{afMax}.png')
    setsProton = buildSets("../out_aseq/Proton*.ASEQ",cov,afMin,afMax)
    plotUpset(setsProton,f'./images/Proton_{cov}_{afMin}_{afMax}.png')

    setTotal = {}
    setTotal["Carbon"] = unionSubset(setsCarbon)
    setTotal["Proton"] = unionSubset(setsProton)
    setTotal["X-ray"] = unionSubset(setsXray)
    setTotal["mock"] = unionSubset(setsMock)

    plotUpset(setTotal, f'./images/allConditions_{cov}_{afMin}_{afMax}.png')

    setTotal["Carbon"] = unionSubsetCarbon(setsCarbon)
    plotUpset(setTotal, f'./images/allConditions_CarbonTOT2Removed_{cov}_{afMin}_{afMax}.png')