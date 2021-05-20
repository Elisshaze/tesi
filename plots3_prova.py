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

afMin = 0.1
afMax = 0.9

out_pacbam = "/media/elisa/backup/out_pacbam/"
out_aseq = "/media/elisa/backup/out_aseq/"
out_aseq_recalib = "/media/elisa/backup/out_aseq_recalib/"

f1= f'{out_aseq}mock_POL_1_S18.final.PILEUP.ASEQ'
f2= f'{out_pacbam}mock_POL_1_S18.final.snps'

d1= pd.read_csv(f1, sep="\t")
d2= pd.read_csv(f2, sep="\t")

df_renamed = d1.rename(columns = {'dbsnp' : 'rsid'})
print(df_renamed)
inter = pd.merge(df_renamed, d2, on="rsid", how="inner")		
filt=((inter['af_x'] >= afMin) & (inter['af_x'] <= afMax)) | ((inter['af_y'] >= afMin) & (inter['af_y'] <= afMax))
inter = inter[filt]

df1=inter[['af_x', 'af_y']]
fig=plt.figure()
sns.kdeplot(data = df1, bw_adjust=.25)
plt.show()
plt.close()

for index in inter.index:
	afx=inter.loc[index, "af_x"]
	afy=inter.loc[index, "af_y"]
	if ((afx >= afMin and afx <= afMax) and (afy < afMin or afy > afMax)):
		inter.loc[index, "color"] = 0
	elif ((afy >= afMin and afy <= afMax) and (afx < afMin or afx > afMax)):
		inter.loc[index, "color"] = 1
	else :
		inter.loc[index, "color"] = 2

fig = plt.figure()
sns.scatterplot(data = inter,  x="af_x", y="af_y", hue="color", palette="hsv", s=40, legend=False)
plt.savefig("why1.png")