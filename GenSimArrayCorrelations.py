import pandas as pd
import numpy as np
import SimulationAlgorithm
import networkx as nx
import time
def basePath():
    return '/slade/home/ss1453/Projects/Other/BackgroundT1DGRS/'
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
    # print(listEither)
    adjacencyMatrix=np.zeros((len(listEither),len(listEither)),float)
    for i,row in tabCorrs.iterrows():
        snp1,snp2=row['SNP1'],row['SNP2']
        if snp1 in listEither and snp2 in listEither:
            idx1,idx2=listEither.index(snp1),listEither.index(snp2)
            adjacencyMatrix[idx1,idx2]=row['Correlation']
            adjacencyMatrix[idx2,idx1]=row['Correlation']
    return adjacencyMatrix,listEither
def returnTableScoresWithCorrelations(size1,dataType,grsType,population,deviationSNPs,tolerance):
    adjacencyMatrix,corrSNPlist=returnAdjacencyMatrix(dataType,grsType,population)
    G=nx.Graph(adjacencyMatrix)
    listConnected=list(nx.connected_components(G)) 
    # print(corrSNPlist)
    # print('Number of groups of SNPs',len(listConnected))
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
            # print('Entering simulation with number of SNPs:',len(scoresList),corrsAll)
            # print(corrsAll)
            # break
            scoreListNew,success=SimulationAlgorithm.returnUpdatedScores(scoresList,size1,corrsAll,tolerance)
            if not success:
                break
            # print('Exiting simulation')
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


deviationSNPsEUR={'rs17843689':[0.803,0.062,0.135]}
deviationSNPsAFR={'rs17843689':[0.69,0.097,0.213]}
deviationSNPsSAS={'rs17843689':[0.963,0.02,0.016]}
deviationSNPsAMR={'rs17843689':[0.885,0.023,0.092]}
deviationSNPsEAS={'rs17843689':[0.901,0.044,0.056]}

deviationSuperPops={'EUR':deviationSNPsEUR,'AFR':deviationSNPsAFR,'SAS':deviationSNPsSAS,
                    'AMR':deviationSNPsAMR,'EAS':deviationSNPsEAS}
tolerance=0.1
maxTime=60
size1=1000
popsWanted=['EUR','AFR','SAS','AMR','EAS']
failedWith=[]
# for population in popsWanted:
#     # population='AFR'
    
#     dataType,grsType='InputDataOther/1000G/','GRS2' #### 
#     deviationSNPs=deviationSuperPops[population]
#     tableFreq=pd.read_csv(path0(dataType,grsType,population)+'tabFreqs'+grsType+population+'.csv',index_col=0)
    
#     tableScores,success=returnTableScoresWithCorrelations(size1,dataType,grsType,population,deviationSNPs,tolerance)
#     time0=time.time()
#     tableFreq=pd.read_csv(path0(dataType,grsType,population)+'tabFreqs'+grsType+population+'.csv',index_col=0)
    
    # while not success:
    #     print('Failed first time, try again')
    #     tableScores,success=returnTableScoresWithCorrelations(size1,dataType,grsType,population,deviationSNPs,tolerance)
    #     elapsedTime=time.time()-time0
    #     if elapsedTime>maxTime:
    #         failedWith.append(population)
    #         break
    # if success:
    #     print('Completed',population)
    # else:
    #     print('Failed',population)
    
    


# size1=200

# listPops=pd.read_csv(path1()+'listPopulations.csv')['Population'].tolist()
# deviantFracsByPop=pd.read_csv(path1()+'populationDeviantSNPfracs.csv',index_col=0)
# for population in listPops:
#     # print(population)
#     # if population=='FIN':
#     dataType,grsType='LDlinkData','GRS2'
#     deviationSNPs0=deviantFracsByPop.loc[population].values
#     deviationSNPs={'rs17843689':[deviationSNPs0[0],deviationSNPs0[1],deviationSNPs0[2]]}
#     tableScores,success=returnTableScoresWithCorrelations(size1,dataType,grsType,population,deviationSNPs)
#     time0=time.time()
#     tableFreq=pd.read_csv(path0(dataType,grsType,population)+'tabFreqs'+grsType+population+'.csv',index_col=0)
    
#     while not success:
#         tableScores,success=returnTableScoresWithCorrelations(size1,dataType,grsType,population,deviationSNPs)
#         elapsedTime=time.time()-time0
#         if elapsedTime>maxTime:
#             failedWith.append(population)
#             break
#     if success:
#         print('Completed',population)
#     else:
#         print('Failed',population)
#     # print('Completed',population)





deviationSNPs={}
## for UKBB
size1=1000
listPops=['T1D']  ##'AFR','SAS','EAS','AFR','EUR',
tolerance=0.01
for population in listPops:
    print(population)
    dataType,grsType='UKBBdata','GRS2'
    deviationSNPs={}
    tableScores,success=returnTableScoresWithCorrelations(size1,dataType,grsType,population,deviationSNPs,tolerance)
    

















