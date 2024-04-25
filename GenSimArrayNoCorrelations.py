import pandas as pd
import numpy as np
import SimulationAlgorithm
def basePath():
    return '/slade/home/ss1453/Projects/Other/BackgroundT1DGRS/'
def path0(dataType,grsType,population):
    return basePath()+'Data/'+dataType+'/'+grsType+'/'+population+'/Tables/'
def path1():
    return basePath()+'Data/OtherData/'
def pathOut(dataType,grsType):
    return basePath()+'SimData/'+dataType+'/'+grsType+'/'
def saveSimSNParray(dataType,grsType,population):
    tableFreq=pd.read_csv(path0(dataType,grsType,population)+'tabFreqs'+grsType+population+'.csv',index_col=0)
    size1=1000
    array1=np.zeros((size1,tableFreq.shape[0]))
    iter0=0
    for snpName,row in tableFreq.iterrows():
        scores=SimulationAlgorithm.returnSNPvals(size1,row['Freq2'])
        array1[:,iter0]=scores
        iter0+=1
    tableScores=pd.DataFrame(columns=list(tableFreq.index),data=array1)
    tableScores.to_csv(pathOut(dataType,grsType)+'SimArrNoCorr'+population+'.csv')
    return tableScores

popsWanted=['AFR','SAS','AMR','EAS'] ##,'EUR'
for population in popsWanted:
    dataType,grsType='LDlinkData','GRS2'
    tableScores=saveSimSNParray(dataType,grsType,population)
listPops=pd.read_csv(path1()+'listPopulations.csv')['Population'].tolist()
for population in listPops:
    dataType,grsType='LDlinkData','GRS2'
    tableScores=saveSimSNParray(dataType,grsType,population)
    print(population)
    
    
## for UKBB
size1=1000
listPops=['EUR','AFR','SAS','EAS']###,'AFR','SAS',

for population in listPops:
    print(population)
    dataType,grsType='UKBBdata','GRS2'
    tableScores=saveSimSNParray(dataType,grsType,population)
       

    
    
    
    
    
    