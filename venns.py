import pandas as pd 
import glob
import os 
import plotly.graph_objects as go
from itertools import cycle
from supervenn import supervenn
import seaborn as sns
import matplotlib.pyplot as plt 
import numpy as np



notre="cov_filtering/notre/"
ricalib="cov_filtering/ricalib/"
res_path="results_cov/venn_pdfs/"

init= ["Carbon_POL", "Carbon_TOT", "mock_POL", "mock_TOT", "Proton_POL", "Proton_TOT", "X-ray_POL", "X-ray_TOT"]

for i in range(2):
	folder = notre
	if i == 1 :
		folder = ricalib

	for b in init :
		fil10 = glob.glob(folder+b+'*.10filtered')
		fil20 = glob.glob(folder+b+'*.20filtered')

		sets10=[]
		for j in fil10:
			df = pd.read_csv(j, sep='\t')
			s = set(df.dbsnp.ravel())
			sets10.append(s)
		
		sets20=[]
		for j in fil10:
			df = pd.read_csv(j, sep='\t')
			s = set(df.dbsnp.ravel())
			sets20.append(s)
		
		pdf= res_path+b+"_"

		f10=[w[w.rindex("/")+1:] for w in fil10]
		f20=[w[w.rindex("/")+1:] for w in fil20]

		
		supervenn(sets10, f10, widths_minmax_ratio=0.1, min_width_for_annotation=900, col_annotations_area_height=2, rotate_col_annotations=True)

		#plt.show()
		
		plt.savefig(pdf+"10fil"+'.pdf')

		plt.close()

		supervenn(sets20, f20, widths_minmax_ratio=0.1, min_width_for_annotation=900, col_annotations_area_height=2, rotate_col_annotations=True,)
		plt.savefig(pdf+"20fil"+'.pdf')
		
		plt.close()
		
