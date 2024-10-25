import pandas as pd
import os
from fnmatch import fnmatch

def basePath():
    return ''
def path0():
    return basePath()+'CorrelationData/'
def pathOut():
    return basePath()+'Tables/'
def saveTableOrCorrelations(grsName,popWanted,minCorr):
    tableCorrelations=pd.DataFrame(columns=['Chr','SNP1','SNP2','Correlation'])
    listNames=os.listdir(path0())
    iter0=0
    for name in listNames:
        file0=pd.read_csv(path0()+name,sep='\t',index_col=0)
        snpNames=list(file0.index)
        if len(snpNames)>1: 
            values=file0.values
            for i in range(values.shape[0]-1):
                for j in range(i+1,values.shape[0]):
                    corr1=values[i,j]
                    if abs(corr1)>minCorr:
                        tableCorrelations.loc[iter0]=[name[0:-4],snpNames[i],snpNames[j],corr1]
                        iter0+=1
    tableCorrelations.to_csv(pathOut()+'tabCorrelations'+grsName+popWanted+'.csv')
    return tableCorrelations

grsName,minCorr='T2DGRS',0.05 ### alter as desired. 
##### examples of running for the 1000G superpopulations and populations
popWanted='EUR'
saveTableOrCorrelations(grsName,popWanted,minCorr)




