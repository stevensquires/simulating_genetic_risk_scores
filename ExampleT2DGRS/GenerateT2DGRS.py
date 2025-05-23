import pandas as pd
import numpy as np
def returnTableScores1(args):
    tableSNPs=pd.read_csv(args['path1']+'SimArr.csv',index_col=0)
    freqData=pd.read_csv(args['path0']+'tabFrequencies.csv',index_col=0)
    scoresFile=pd.read_excel(args['path0']+'grsScores.xlsx').set_index('RSID')
    allAvailSNPs=list(set(freqData.index.tolist()).intersection(set(scoresFile.index.tolist())))
    print(len(allAvailSNPs))
    scoresFile=scoresFile.loc[allAvailSNPs]
    arraySNPs2=np.zeros((tableSNPs.shape[0],tableSNPs.shape[1]),int)
    arrScores1=np.zeros((tableSNPs.shape[0],tableSNPs.shape[1]),float)
    snpNames,iter0=[],0
    for snpName,scoreRow in scoresFile.iterrows():
        freqTemp=freqData.loc[snpName]
        vals=tableSNPs[snpName].values
        effectAllele= scoreRow['A1']  ### score alleles are A1 for this T2DGRS
        if effectAllele==freqTemp['A1']:
            vals2=2-vals
        elif effectAllele==freqTemp['A2']:
            vals2=vals
        else:
            print(snpName) ## check no issue with alleles
        arraySNPs2[:,iter0]=vals2
        arrScores1[:,iter0]=vals2*scoreRow['BETA']
        snpNames.append(snpName)
        iter0+=1
    tableScores=pd.DataFrame(columns=snpNames,data=arrScores1)
    tableSNPs2=pd.DataFrame(columns=snpNames,data=arraySNPs2)
    return tableScores,tableSNPs2,scoresFile

    
def saveAllScores(args):
    tableScores,tableSNPs2,scoresFile=returnTableScores1(args)
    summedScore=np.sum(tableScores.values,axis=1)
    np.save(args['path1']+'FinalScores.npy',summedScore)
    tableFinalScores=pd.DataFrame(columns=['FinalScore'],data=summedScore)
    tableFinalScores.to_csv(args['path1']+'FinalScores.csv')
    tableScores.to_csv(args['path1']+'ScoresPerSNP.csv')
    return summedScore





