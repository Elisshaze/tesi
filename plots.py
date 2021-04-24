import pandas as pd 
import os 
import plotly.graph_objects as go

Afilt=[]
Bfilt=[]
Cfilt=[]
Dfilt=[]

#file che contiene il numero di righe in merge e intersezioni avvenute e il numero di righe del dataframe

stats = open("stats.txt", "a")


colnames = ['chr','pos','dbsnp','MAF','ref','alt','A','C','G','T','RD','Ars','Crs','Grs','Trs','af','cov']

#######creazione dei dataframe per pileup non ricalibrati####################################
directory="local1/"
for filename in os.listdir(directory):
    if filename.endswith(".Afiltered"): 
        fil=os.path.join(directory, filename)
        dfA = pd.read_csv(fil, names=colnames, sep='\t')
        
        Afilt.append(dfA) 
    if filename.endswith(".Bfiltered"): 
        fil=os.path.join(directory, filename)
        dfB = pd.read_csv(fil, names=colnames, sep='\t')
        
        Bfilt.append(dfB)
    if filename.endswith(".Cfiltered"): 
        fil=os.path.join(directory, filename)
        dfC = pd.read_csv(fil, names=colnames, sep='\t')
        
        Cfilt.append(dfC)
    if filename.endswith(".Dfiltered"): 
        fil=os.path.join(directory, filename)
        dfD = pd.read_csv(fil, names=colnames, sep='\t')
        
        Dfilt.append(dfD)

###########################################################################
#merge and intersect of all files of the same filtering for recalibrated
##########################################################################

boh=[Afilt, Bfilt, Cfilt, Dfilt]

names=["Afilt", "Bfilt", "Cfilt", "Dfilt"]

j=0
for filt in boh :
    filt_merge=filt[0]
    for df in filt[1:]:
        filt_merge = pd.merge(filt_merge, df, on='dbsnp', how='inner')
   
    filt_merge=pd.concat(filt).drop_duplicates().reset_index(drop=True)

    stats.write("merge"+'\t'+names[j]+'\t'+str(len(filt_merge.index))+'\n')
    filename="csvs/"+names[j]+".csv"

    open(filename, "w")
    filt_merge.to_csv(filename, sep='\t')
    j+=1

#######creazione dei dataframe per pileup ricalibrati####################################

ARfilt=[]
BRfilt=[]
CRfilt=[]
DRfilt=[]

directory="local/"
for filename in os.listdir(directory):
    if filename.endswith(".Afiltered"): 
        fil=os.path.join(directory, filename)
        dfA = pd.read_csv(fil, names=colnames, sep='\t')
        
        ARfilt.append(dfA) 
    if filename.endswith(".Bfiltered"): 
        fil=os.path.join(directory, filename)
        dfB = pd.read_csv(fil, names=colnames, sep='\t')
        
        BRfilt.append(dfB)
    if filename.endswith(".Cfiltered"): 
        fil=os.path.join(directory, filename)
        dfC = pd.read_csv(fil, names=colnames, sep='\t')
        
        CRfilt.append(dfC)
    if filename.endswith(".Dfiltered"): 
        fil=os.path.join(directory, filename)
        dfD = pd.read_csv(fil, names=colnames, sep='\t')
        
        DRfilt.append(dfD)

#############################################################################
#merge and intersect of all files of the same filtering for non recalibrated
#############################################################################


boh=[ARfilt, BRfilt, CRfilt, DRfilt]
names=["ARfilt", "BRfilt", "CRfilt", "DfRfilt"]

j=0
for filt in boh :
    filt_merge=filt[0]
    for df in filt[1:]:
        filt_merge = pd.merge(filt_merge, df, on='dbsnp', how='inner')
   
    filt_merge=pd.concat(filt).drop_duplicates().reset_index(drop=True)

    stats.write("merge"+'\t'+names[j]+'\t'+str(len(filt_merge.index))+'\n')
    filename="csvs/"+names[j]+".csv"

    open(filename, "w")
    filt_merge.to_csv(filename, sep='\t')
    j+=1


######################################################################################
#SCATTER PLOTS BETWEEN RECALIBRATED AND NOT RECALIBRATED MOCKS FOR DIFFERENT REPLICATE
#####################################################################################

#name of files, only files with cov >20 and 0.2 < af < 0.8 
mock_TOT_1_R='local/mock_TOT_1_S15.final.groups.dedup.splitted.realigned.recalibrated.PILEUP.ASEQ.Dfiltered'
mock_TOT_3_R='local/mock_TOT_3_S14.final.groups.dedup.splitted.realigned.recalibrated.PILEUP.ASEQ.Dfiltered'
mock_TOT_6_R='local/mock_TOT_6_S13.final.groups.dedup.splitted.realigned.recalibrated.PILEUP.ASEQ.Dfiltered'

mock_TOT_1='local1/mock_TOT_1_S15.final.PILEUP.ASEQ.Dfiltered'
mock_TOT_3='local1/mock_TOT_3_S14.final.PILEUP.ASEQ.Dfiltered'
mock_TOT_6='local1/mock_TOT_6_S13.final.PILEUP.ASEQ.Dfiltered'

#creation of dataframes
df1TR = pd.read_csv(mock_TOT_1_R, names=colnames, sep='\t')
df3TR = pd.read_csv(mock_TOT_3_R, names=colnames, sep='\t')
df6TR = pd.read_csv(mock_TOT_6_R, names=colnames, sep='\t')

df1T = pd.read_csv(mock_TOT_1, names=colnames, sep='\t')
df3T = pd.read_csv(mock_TOT_3, names=colnames, sep='\t')
df6T = pd.read_csv(mock_TOT_6, names=colnames, sep='\t')

#lists containing dataframes
rec_array=[df1TR, df3TR, df6TR]
array=[df1T, df3T, df6T]

#intersection and union of data frames, common SNPs for mocks
#ALL RECALIBRATED
mocks_rec_merge = rec_array[0]
for df in rec_array[1:]:       
    mocks_rec_merge = pd.merge(mocks_rec_merge, df, on='dbsnp', how='inner')

open("mocks_rec_intersect.csv","w")
mocks_rec_merge.to_csv("mocks_rec_intersect.csv", sep='\t')
stats.write("inter"+'\t'+"all_rec"+'\t'+str(len(mocks_rec_merge))+'\n')

mocks_rec_merge=pd.concat(rec_array).drop_duplicates().reset_index(drop=True)
open("mocks_rec_merge.csv","w")
mocks_rec_merge.to_csv("mocks_rec_merge.csv", sep='\t')
stats.write("merge"+'\t'+"all_rec"+'\t'+str(len(mocks_rec_merge))+'\n')

#ALL NOT RECALIBRATED
mocks_merge = array[0]
for df in array[1:]:       
    mocks_merge = pd.merge(mocks_merge, df, on='dbsnp', how='inner')

open("csvs/mocks_intersect.csv","w")
mocks_merge.to_csv("csvs/mocks_intersect.csv", sep='\t')
stats.write("inter"+'\t'+"all_not"+'\t'+str(len(mocks_merge))+'\n')

mocks_merge=pd.concat(array).drop_duplicates().reset_index(drop=True)
open("csvs/mocks_merge.csv","w")
mocks_merge.to_csv("csvs/mocks_merge.csv", sep='\t')
stats.write("merge"+'\t'+"all_not"+'\t'+str(len(mocks_merge))+'\n')


####################################################################################
#Intersection between recalibrated pileups
TR1_R3=pd.merge(df1TR, df3TR, on="dbsnp", how='inner')
stats.write("inter"+'\t'+"TR1_R3"+'\t'+str(len(TR1_R3))+'\n')
open("csvs/TR1_R3.csv","w")
TR1_R3.to_csv("csvs/TR1_R3.csv", sep='\t')

TR1_R6=pd.merge(df1TR, df6TR, on="dbsnp", how='inner')
stats.write("inter"+'\t'+"TR1_R6"+'\t'+str(len(TR1_R6))+'\n')
open("csvs/TR1_R6.csv","w")
TR1_R6.to_csv("csvs/TR1_R6.csv", sep='\t')

TR3_R6=pd.merge(df6TR, df3TR, on="dbsnp", how='inner')
stats.write("inter"+'\t'+"TR3_R6"+'\t'+str(len(TR3_R6))+'\n')
open("csvs/TR3_R6.csv","w")
TR3_R6.to_csv("csvs/TR3_R6.csv", sep='\t')

####################################################################################
#intersection between non recalibrated pileups
T1_3=pd.merge(df1T, df3T, on="dbsnp", how='inner')
stats.write("inter"+'\t'+"T1_3"+'\t'+str(len(T1_3))+'\n')
open("csvs/T1_3.csv","w")
T1_3.to_csv("csvs/T1_3.csv", sep='\t')

T1_6=pd.merge(df1T, df6T, on="dbsnp", how='inner')
stats.write("inter"+'\t'+"T1_6"+'\t'+str(len(T1_6))+'\n')
open("csvs/T1_6.csv","w")
T1_6.to_csv("csvs/T1_6.csv", sep='\t')

T3_6=pd.merge(df6T, df3T, on="dbsnp", how='inner')
stats.write("inter"+'\t'+"T3_6"+'\t'+str(len(T3_6))+'\n')
open("csvs/T3_6.csv","w")
T3_6.to_csv("csvs/T3_6.csv", sep='\t')

#####################################################################################
#intersection between recalibrated and non recalibrated pileups
TR1_1=pd.merge(df1TR, df1T, on="dbsnp", how='inner')
stats.write("inter"+'\t'+"TR1_1"+'\t'+str(len(TR1_1))+'\n')
open("csvs/TR1_1.csv","w")
TR1_1.to_csv("csvs/TR1_1.csv", sep='\t')

TR3_3=pd.merge(df3TR, df3T, on="dbsnp", how='inner')
stats.write("inter"+'\t'+"TR3_3"+'\t'+str(len(TR3_3))+'\n')
open("csvs/TR3_3.csv","w")
TR3_3.to_csv("csvs/TR3_3.csv", sep='\t')

TR6_6=pd.merge(df6TR, df6T, on="dbsnp", how='inner')
stats.write("inter"+'\t'+"TR6_6"+'\t'+str(len(TR6_6))+'\n')
open("csvs/TR6_6.csv","w")
TR6_6.to_csv("csvs/TR6_6.csv", sep='\t')

#########################################################################################
#CREAZIONE DEI PLOT
###########################################################################################

plot_array=[TR1_R3, TR1_R6, TR3_R6, T1_3, T1_6, T3_6, TR1_1, TR3_3, TR6_6]
plot_names=["TR1_R3", "TR1_R6", "TR3_R6", "T1_3", "T1_6", "T3_6", "TR1_1", "TR3_3", "TR6_6"]

for i in range(len(plot_array)):
    plot="scatter_"+plot_names[i]

    line_x = [0, 1]
    line_y = [0, 1]
    fig1 = go.Figure()
    # Scatter plot of recalibrated of mock_TOT_1 and mock_TOT_3
    fig1.add_trace(go.Scatter(x=plot_array[i]['af_x'], y=plot_array[i]['af_y'], text=plot_array[i]['dbsnp'],
                        mode='markers',
                        name='SNPS'))

    fig1.add_trace(go.Scatter(x=line_x, y=line_y,
                        mode='lines',
                        name='bisector',
                        line = dict(dash = 'dash')))

    fig1.write_html("htmls/"+plot+".html")
    fig1.write_image("pdfs/"+plot+".pdf")
