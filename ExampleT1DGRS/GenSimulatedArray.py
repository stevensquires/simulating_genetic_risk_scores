import pandas as pd
import numpy as np
import SimulationAlgorithm
import networkx as nx
import time

def returnSimDict(args,deviationSNPs):
    tableFreq=pd.read_csv(args['path0']+'tabFrequencies.csv',index_col=0)
    dictScores={}
    for snpName,row in tableFreq.iterrows():
        if snpName in deviationSNPs:
            scores=SimulationAlgorithm.returnSNPvalsIfDeviate(args['number'],deviationSNPs[snpName])
        else:
            scores=SimulationAlgorithm.returnSNPvals(args['number'],row['Freq2'])
        dictScores[snpName]=scores
    return dictScores
def returnAdjacencyMatrix(args):
    tableFreq=pd.read_csv(args['path0']+'tabFrequencies.csv',index_col=0)
    tabCorrs=pd.read_csv(args['path0']+'tabCorrelations.csv',index_col=0)
    list0,list1=set(list(tabCorrs['SNP1'])),set(list(tabCorrs['SNP2']))
    listEither0,listEither=list(list0.union(list1)),[]
    for snpName in listEither0:
        row=tableFreq.loc[snpName]
        freq1,freq2=row['Freq1'],row['Freq2']
        minFreq=np.minimum(freq1,freq2)
        if minFreq>args['minFreq']:
            listEither.append(snpName)
    adjacencyMatrix=np.zeros((len(listEither),len(listEither)),float)
    for i,row in tabCorrs.iterrows():
        snp1,snp2=row['SNP1'],row['SNP2']
        if snp1 in listEither and snp2 in listEither:
            idx1,idx2=listEither.index(snp1),listEither.index(snp2)
            adjacencyMatrix[idx1,idx2]=row['Correlation']
            adjacencyMatrix[idx2,idx1]=row['Correlation']
    return adjacencyMatrix,listEither
def returnTableScoresWithCorrelations(args,deviationSNPs):
    adjacencyMatrix,corrSNPlist=returnAdjacencyMatrix(args)
    G=nx.Graph(adjacencyMatrix)
    listConnected=list(nx.connected_components(G)) 
    # print(listConnected)
    dictScores=returnSimDict(args,deviationSNPs)
    anyFailures=False
    for i,dictConnect1 in enumerate(listConnected):
        listConnect1=list(dictConnect1)
        if len(listConnect1)>1:
            #print(listConnect1)
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
            scoreListNew,success=SimulationAlgorithm.returnUpdatedScores(scoresList,corrsAll,args)
            if not success:
                print('Failed with',listConnect1)
                anyFailures=True
            for i,snpName in enumerate(snpNames):
                dictScores[snpName]=scoreListNew[i]
    snpNamesAll,arrayAll=[],np.zeros((args['number'],len(dictScores)))
    iter0=0
    for snpName,arr1 in dictScores.items():
        arrayAll[:,iter0]=arr1
        snpNamesAll.append(snpName)
        iter0+=1
    tableScores=pd.DataFrame(columns=snpNamesAll,data=arrayAll)
    if success:
        tableScores.to_csv(args['path1']+'SimArr.csv')
    return tableScores,anyFailures
def returnTableScoresWithoutCorrelations(args,deviationSNPs):
    dictScores=returnSimDict(args,deviationSNPs)
    snpNamesAll,arrayAll=[],np.zeros((args['number'],len(dictScores)))
    iter0=0
    for snpName,arr1 in dictScores.items():
        arrayAll[:,iter0]=arr1
        snpNamesAll.append(snpName)
        iter0+=1
    tableScores=pd.DataFrame(columns=snpNamesAll,data=arrayAll)
    tableScores.to_csv(args['path1']+'SimArrNoCorr.csv')
    return tableScores

def genSimArr(args,deviationSNPs):
    tableScores,success=returnTableScoresWithCorrelations(args,deviationSNPs)
    return tableScores,success








