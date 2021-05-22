import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

lst = [[3671023, 0.200000, 0.333333], [4492071, 0.176471, 0.333333],
       [4492302, 0.222222, 0.285714], [4525905, 0.298246, 0.234043],
	   [4520905, 0.003334, 0.234043], [4520905, 0.400098, 0.000221], 
	   [4520905, 0.001134, 0.714043], [4520905, 0.559008, 0.010221]
	   ]
df = pd.DataFrame(lst, columns =['pos', 'af_x', 'af_y'])

afMin=0.1
afMax=0.9

df['color']=np.nan
for index in df.index:
	afx=df.loc[index, "af_x"]
	afy=df.loc[index, "af_y"]
	if ((afx >= afMin and afx <= afMax) and (afy < afMin or afy > afMax)):
		df.loc[index, "color"] = 0
	elif ((afy >= afMin and afy <= afMax) and (afx < afMin or afx > afMax)):
		df.loc[index, "color"] = 1
	elif ((afy >= afMin and afy <= afMax) and (afx >= afMin or afx <= afMax)):
		df.loc[index, "color"] = 2

sns.scatterplot(data = df,  x="af_x", y="af_y", hue="color", palette = "hsv", s=40, legend=False)

plt.savefig("stack_why_hsv.png")

		#violet		#green		#orange
colors = ['#747FE3', '#8EE35D', '#E37346']
sns.set_palette(sns.color_palette(colors))

sns.scatterplot(data = df,  x="af_x", y="af_y", hue="color", s=40, legend=False)
plt.savefig("stack_why_personal.png")