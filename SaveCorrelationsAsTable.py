import pandas as pd
import os
from fnmatch import fnmatch

def basePath():
    return '/slade/home/ss1453/Projects/Other/BackgroundT1DGRS/'
def path0(grsName,popWanted):
    return basePath()+'Data/LDlinkData/'+grsName+'/'+popWanted+'/Raw/'
def path1():
    return basePath()+'Data/OtherData/'
def pathOut(grsName,popWanted):
    return basePath()+'Data/LDlinkData/'+grsName+'/'+popWanted+'/Tables/'
def returnCorreNames(grsName,popWanted):
    list0=os.listdir(path0(grsName,popWanted))
    listChrs=[]
    for name in list0:
        if fnmatch(name,'*Correlations*'):
            listChrs.append(name.split('Correlations')[0])
    return list(set(listChrs))
def saveTableOrCorrelations(grsName,popWanted,minCorr):
    tableCorrelations=pd.DataFrame(columns=['Chr','SNP1','SNP2','Correlation'])
    listNames=returnCorreNames(grsName,popWanted)
    iter0=0
    for chrName in listNames:
        file0=pd.read_csv(path0(grsName,popWanted)+chrName+'Correlations.txt',sep='\t',index_col=0)
        snpNames=list(file0.index)
        if len(snpNames)>1:
            values=file0.values
            for i in range(values.shape[0]-1):
                for j in range(i+1,values.shape[0]):
                    corr1=values[i,j]
                    if abs(corr1)>minCorr:
                        tableCorrelations.loc[iter0]=[chrName,snpNames[i],snpNames[j],corr1]
                        iter0+=1
    tableCorrelations.to_csv(pathOut(grsName,popWanted)+'tabCorrelations'+grsName+popWanted+'.csv')
    return tableCorrelations

popsWanted=['AFR','SAS','AMR','EAS']
for popWanted in popsWanted:
    grsName,minCorr='GRS2',0.05
    tableCorrelations=saveTableOrCorrelations(grsName,popWanted,minCorr)
listPops=pd.read_csv(path1()+'listPopulations.csv')['Population'].tolist()
for popWanted in listPops:
    grsName,minCorr='GRS2',0.05
    tableCorrelations=saveTableOrCorrelations(grsName,popWanted,minCorr)






