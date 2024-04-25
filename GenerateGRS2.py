import pandas as pd
import numpy as np
import pickle
def basePath():
    return '/slade/home/ss1453/Projects/Other/BackgroundT1DGRS/'
def path0(dataType,grsType):
    return basePath()+'SimData/'+dataType+'/'+grsType+'/'
def path1():
    return basePath()+'Data/LDlinkData/'
def path2(dataType,grsType,population):
    return basePath()+'Data/'+dataType+'/'+grsType+'/'+population+'/Tables/'
def path3():
    return basePath()+'Data/OtherData/'
def pathOut():
    return basePath()+'Outputs/Outputs2/'

def returnTableScores1(dataType,grsType,population,simOrNot):
    tableSNPs=pd.read_csv(path0(dataType,grsType)+'SimArr'+simOrNot+population+'.csv',index_col=0)
    freqData=pd.read_csv(path2(dataType,grsType,population)+'tabFreqs'+grsType+population+'.csv',index_col=0)
    if dataType=='LDlinkData':
        scoresFile=pd.read_excel(path1()+'T1DGRS67_1000G_hg19_FinalFor1000GSimulation.xlsx').set_index('RSID')
    elif dataType=='UKBBdata':
        scoresFile=pd.read_excel(path1()+'T1DGRS67_1000G_hg19_FinalFor1000GSimulationUKBB.xlsx').set_index('RSID')
    arraySNPs2=np.zeros((tableSNPs.shape[0],tableSNPs.shape[1]),int)
    arrScores1=np.zeros((tableSNPs.shape[0],tableSNPs.shape[1]),float)
    snpNames,iter0=[],0
    for snpName,scoreRow in scoresFile.iterrows():
        freqTemp=freqData.loc[snpName]
        vals=tableSNPs[snpName].values
        effectAllele= scoreRow['SCORE_ALLELE']
        if effectAllele==freqTemp['A1']:
            vals2=2-vals
        elif effectAllele==freqTemp['A2']:
            vals2=vals
        else:
            print(snpName) ## check no issue with alleles
        arraySNPs2[:,iter0]=vals2
        arrScores1[:,iter0]=vals2*scoreRow['SCORE']
        snpNames.append(snpName)
        iter0+=1
    tableScores=pd.DataFrame(columns=snpNames,data=arrScores1)
    tableSNPs2=pd.DataFrame(columns=snpNames,data=arraySNPs2)
    return tableScores,tableSNPs2,scoresFile
def returnDQscores(tableSNPs2,scoresFile):
    intScores=pd.read_csv(path1()+'T1D67_int_scores.txt',sep='\t')
    ranking=pd.read_csv(path1()+'HLADQ_USAEuropean_Klitz.txt',sep='\t')
    rank1=list(ranking['DQ'])
    listDQ,listDQsorted,dqScores=[],[],[]
    arrSNPs=tableSNPs2.values
    for i in range(arrSNPs.shape[0]):
        listTemp=[]
        iter0=0
        for snpName,scoreRow in scoresFile.iterrows():
            val=arrSNPs[i,iter0]
            dqName=scoreRow['COMPONENT'].split('-')[1]
            if val==1:
                listTemp.append(dqName)
            elif val==2:
                listTemp.append(dqName)
                listTemp.append(dqName)
            iter0+=1
            if iter0==14:
                break
        listDQ.append(listTemp)
        temp2 = [val1 for x in rank1 for val1 in listTemp if val1 == x]
        listDQsorted.append(temp2)
        dqScore=0
        for i,row in intScores.iterrows():
            if len(temp2)>1:
                if row['ALLELE1']==temp2[0] and row['ALLELE2']==temp2[1]:
                    dqScore=row['BETA']
                    break
        dqScores.append(dqScore)
    return np.array(dqScores),listDQ,listDQsorted
    
def saveAllScores(dataType,grsType,population,corrOrNot):
    tableScores,tableSNPs2,scoresFile=returnTableScores1(dataType,grsType,population,corrOrNot)
    dqScores,listDQ,listDQsorted=returnDQscores(tableSNPs2,scoresFile)
    summedScore=np.sum(tableScores.values,axis=1)
    finalScores=summedScore+dqScores
    np.save(pathOut()+'FinalScores'+dataType+grsType+population+corrOrNot+'.npy',finalScores)
    
    arrayAll=np.concatenate([finalScores.reshape(-1,1),summedScore.reshape(-1,1),dqScores.reshape(-1,1)],axis=1)
    
    tableFinalScores=pd.DataFrame(columns=['FinalScore','SumScore','InteractionScore'],data=arrayAll)
    tableFinalScores.to_csv(pathOut()+'FinalScores'+dataType+grsType+population+corrOrNot+'.csv')
    pickle.dump(listDQsorted,open(pathOut()+'DQlistSorted'+dataType+grsType+population+corrOrNot+'.pkl','wb'))
    tableScores.to_csv(pathOut()+'ScoresPerSNP'+dataType+grsType+population+corrOrNot+'.csv')
    return tableFinalScores


# popsWanted=['EUR','AFR','SAS','AMR','EAS'] ##,'EUR'
# corrOrNot=['Corr']#,'NoCorr']##'NoCorr' ### 
# for population in popsWanted:
#     for corr in corrOrNot:
#         dataType,grsType='LDlinkData','GRS2'
#         tableFinalScores=saveAllScores(dataType,grsType,population,corr)
        
    
# listPops=pd.read_csv(path3()+'listPopulations.csv')['Population'].tolist()
# corrOrNot=['Corr']#,'Corr']###,'NoCorr']##'NoCorr' ### 
# for population in listPops:
#     for corr in corrOrNot:
#         dataType,grsType='LDlinkData','GRS2'
#         tableFinalScores=saveAllScores(dataType,grsType,population,corr)
#     print(population)
      

    
listPops=['T1D']#,'AFR','SAS','EAS','EUR'] 
corrOrNot=['Corr']#,'Corr']###,'NoCorr']##'NoCorr' ### 
for population in listPops:
    for corr in corrOrNot:
        dataType,grsType='UKBBdata','GRS2'
        tableFinalScores=saveAllScores(dataType,grsType,population,corr)
    print(population)



