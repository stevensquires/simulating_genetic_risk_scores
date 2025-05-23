import pandas as pd
import numpy as np
import SimulationAlgorithm
import networkx as nx
import time
def basePath():
    return ''
def path0(dataType,grsType,population):
    return basePath()+'Data/'+dataType+'/'+grsType+'/'+population+'/Tables/'
def path1():
    return basePath()+'Data/OtherData/'
def pathOut(dataType,grsType):
    return basePath()+'SimData/'+dataType+'/'+grsType+'/'
def returnSimDict(dataType,grsType,population,deviationSNPs,size1):
    tableFreq=pd.read_csv(path0(dataType,grsType,population)+'tabFreqs'+grsType+population+'.csv',index_col=0)
    dictScores={}
    for snpName,row in tableFreq.iterrows():
        if snpName in deviationSNPs:
            scores=SimulationAlgorithm.returnSNPvalsIfDeviate(size1,deviationSNPs[snpName])
        else:
            scores=SimulationAlgorithm.returnSNPvals(size1,row['Freq2'])
        dictScores[snpName]=scores
    return dictScores
def returnAdjacencyMatrix(dataType,grsType,population):
    tableFreq=pd.read_csv(path0(dataType,grsType,population)+'tabFreqs'+grsType+population+'.csv',index_col=0)
    tabCorrs=pd.read_csv(path0(dataType,grsType,population)+'tabCorrelations'+grsType+population+'.csv',index_col=0)
    list0,list1=set(list(tabCorrs['SNP1'])),set(list(tabCorrs['SNP2']))
    listEither0,listEither=list(list0.union(list1)),[]
    for snpName in listEither0:
        row=tableFreq.loc[snpName]
        freq1,freq2=row['Freq1'],row['Freq2']
        minFreq=np.minimum(freq1,freq2)
        if minFreq>0.01:
            listEither.append(snpName)
    adjacencyMatrix=np.zeros((len(listEither),len(listEither)),float)
    for i,row in tabCorrs.iterrows():
        snp1,snp2=row['SNP1'],row['SNP2']
        if snp1 in listEither and snp2 in listEither:
            idx1,idx2=listEither.index(snp1),listEither.index(snp2)
            adjacencyMatrix[idx1,idx2]=row['Correlation']
            adjacencyMatrix[idx2,idx1]=row['Correlation']
    return adjacencyMatrix,listEither
def returnTableScoresWithCorrelations(size1,dataType,grsType,population,deviationSNPs,tolerance,iterationLimit):
    adjacencyMatrix,corrSNPlist=returnAdjacencyMatrix(dataType,grsType,population)
    G=nx.Graph(adjacencyMatrix)
    listConnected=list(nx.connected_components(G)) 
    # print(listConnected)
    dictScores=returnSimDict(dataType,grsType,population,deviationSNPs,size1)
    for i,dictConnect1 in enumerate(listConnected):
        listConnect1=list(dictConnect1)
        if len(listConnect1)>1:
            print(listConnect1)
            corrsAll,scoresList,snpNames={},[],[]
            for idx in listConnect1:
                snpNames.append(corrSNPlist[idx])
                scoresList.append(dictScores[corrSNPlist[idx]])
            for i in range(len(listConnect1)-1):
                idx0=listConnect1[i]
                for j in range(i+1,len(listConnect1)):
                    idx1=listConnect1[j]
                    corrTemp=adjacencyMatrix[idx0,idx1]
                    if corrTemp!=0:
                        corrsAll[i,j]=adjacencyMatrix[idx0,idx1]
            scoreListNew,success=SimulationAlgorithm.returnUpdatedScores(scoresList,size1,corrsAll,tolerance,iterationLimit)
            if not success:
                print('Failed with',listConnect1)
            for i,snpName in enumerate(snpNames):
                dictScores[snpName]=scoreListNew[i]
    snpNamesAll,arrayAll=[],np.zeros((size1,len(dictScores)))
    iter0=0
    for snpName,arr1 in dictScores.items():
        arrayAll[:,iter0]=arr1
        snpNamesAll.append(snpName)
        iter0+=1
    tableScores=pd.DataFrame(columns=snpNamesAll,data=arrayAll)
    if success:
        tableScores.to_csv(pathOut(dataType,grsType)+'SimArrCorr'+population+'.csv')
    return tableScores,success
def returnTableScoresWithoutCorrelations(dataType,grsType,population,deviationSNPs,size1):
    dictScores=returnSimDict(dataType,grsType,population,deviationSNPs,size1)
    snpNamesAll,arrayAll=[],np.zeros((size1,len(dictScores)))
    iter0=0
    for snpName,arr1 in dictScores.items():
        arrayAll[:,iter0]=arr1
        snpNamesAll.append(snpName)
        iter0+=1
    tableScores=pd.DataFrame(columns=snpNamesAll,data=arrayAll)
    tableScores.to_csv(pathOut(dataType,grsType)+'SimArrNoCorr'+population+'.csv')
    return tableScores


deviationSNPs={'rs17843689':[0.803,0.062,0.135]} ### example of a SNP out of HWE that needs to be generated separately

dataType,grsType='LDlinkData','GRS2'  #### can be changed as desired - alters the directories
population='EUR' ### relevant for 1000G data, may not be relevant for others - alters the directories
size1=1000 ### number of simulated samples to be generated

wantSNParrayWithoutCorrelations=True
if wantSNParrayWithoutCorrelations:
    tableScoresNoCorrelations=returnTableScoresWithoutCorrelations(dataType,grsType,population,deviationSNPs,size1)
tolerance=0.01
iterationLimit=10000
wantSNParrayWithCorrelations=True
if wantSNParrayWithCorrelations:
    tableScoresWithCorrelations=returnTableScoresWithCorrelations(size1,dataType,grsType,population,deviationSNPs,tolerance,iterationLimit)












