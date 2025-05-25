import pandas as pd
import numpy as np
import pickle

def returnTableScores1(args):
    tableSNPs=pd.read_csv(args['path1']+'SimArr.csv',index_col=0)
    freqData=pd.read_csv(args['path0']+'tabFrequencies.csv',index_col=0)
    scoresFile=pd.read_excel(args['path0']+'grsScores.xlsx').set_index('RSID')
    allAvailSNPs=freqData.index.tolist()
    scoresFile=scoresFile.loc[allAvailSNPs]
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
def returnDQscores(tableSNPs2,scoresFile,args):
    intScores=pd.read_csv(args['path0']+'interaction_scores.txt',sep='\t')
    ranking=pd.read_csv(args['path0']+'rankingFile.txt',sep='\t')
    rank1=list(ranking['DQ'])
    listDQ,listDQsorted,dqScores=[],[],[]
    arrSNPs=tableSNPs2.values
    
    for i in range(arrSNPs.shape[0]):
        listTemp=[]
        iter0=0
        for snpName,scoreRow in scoresFile.iterrows():
            val=arrSNPs[i,iter0]
            dq=scoreRow['COMPONENT']
            if 'DQ' in dq:
                dqName=scoreRow['COMPONENT'].split('-')[1]
                if val==1:
                    listTemp.append(dqName)
                elif val==2:
                    listTemp.append(dqName)
                    listTemp.append(dqName)
            iter0+=1

        listDQ.append(listTemp)
        temp2 = [val1 for x in rank1 for val1 in listTemp if val1 == x]
        listDQsorted.append(temp2)
        dqScore=0
        for i,row in intScores.iterrows():
            if len(temp2)>1:
                if row['ALLELE1']==temp2[0] and row['ALLELE2']==temp2[1]:
                    dqScore=row['BETA']
                    break
                if row['ALLELE2']==temp2[0] and row['ALLELE1']==temp2[1]:
                    dqScore=row['BETA']
                    break
        dqScores.append(dqScore)
    return np.array(dqScores),listDQ,listDQsorted
    
def saveAllScores(args):
    tableScores,tableSNPs2,scoresFile=returnTableScores1(args)
    dqScores,listDQ,listDQsorted=returnDQscores(tableSNPs2,scoresFile,args)
    summedScore=np.sum(tableScores.values,axis=1)
    finalScores=summedScore+dqScores
    np.save(args['path1']+'FinalScores.npy',finalScores)
    
    arrayAll=np.concatenate([finalScores.reshape(-1,1),summedScore.reshape(-1,1),dqScores.reshape(-1,1)],axis=1)
    
    tableFinalScores=pd.DataFrame(columns=['FinalScore','SumScore','InteractionScore'],data=arrayAll)
    tableFinalScores.to_csv(args['path1']+'FinalScores.csv')
    pickle.dump(listDQsorted,open(args['path1']+'DQlistSorted.pkl','wb'))
    tableScores.to_csv(args['path1']+'ScoresPerSNP.csv')
    return tableFinalScores




