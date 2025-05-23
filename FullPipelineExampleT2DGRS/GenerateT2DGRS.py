import pandas as pd
import numpy as np
import pickle
def basePath():
    return ''
def path0():
    return basePath()+'SimData/'
def path2():
    return basePath()+'Tables/'
def pathOut():
    return basePath()+'Outputs/'

def returnTableScores1(dataType,grsType,population,simOrNot):
    tableSNPs=pd.read_csv(path0()+'SimArr'+simOrNot+population+'.csv',index_col=0)
    freqData=pd.read_csv(path2()+'tabFreqs'+grsType+population+'.csv',index_col=0)
    allAvailSNPs=freqData.index.tolist()
    scoresFile=pd.read_excel(basePath()+'T2D_GRS392_SEARCH_hg19.xlsx').set_index('RSID')
    arraySNPs2=np.zeros((tableSNPs.shape[0],tableSNPs.shape[1]),int)
    arrScores1=np.zeros((tableSNPs.shape[0],tableSNPs.shape[1]),float)
    snpNames,iter0=[],0
    noPreds=pd.read_csv('plink.nopred',sep='\s+',header=None)[1].tolist()
    noPreds0=[]
    for name in noPreds:
        temp0=name.split(':')
        chrPos=temp0[0]+':'+temp0[1]
        noPreds0.append(chrPos)
        
    for snpName,scoreRow in scoresFile.iterrows():
        chrPos=scoreRow['CHR:POS']#.split(':')
        print(chrPos)
        if chrPos not in noPreds0:
            if snpName in allAvailSNPs:
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
            else:
                print('Missing',snpName)
        else:
            print(chrPos)
    tableScores=pd.DataFrame(columns=snpNames,data=arrScores1)
    tableSNPs2=pd.DataFrame(columns=snpNames,data=arraySNPs2)
    return tableScores,tableSNPs2,scoresFile

    
def saveAllScores(dataType,grsType,population,corrOrNot):
    tableScores,tableSNPs2,scoresFile=returnTableScores1(dataType,grsType,population,corrOrNot)
    summedScore=np.sum(tableScores.values,axis=1)
    np.save(pathOut()+'FinalScores'+dataType+grsType+population+corrOrNot+'.npy',summedScore)
    return summedScore

if __name__=='__main__':
    corrOrNot=['Corr','NoCorr']##'NoCorr' ###
    population,dataType,grsType='EUR','LDlinkData','T2DGRS'
    tableEURnoCorr=tableFinalScores=saveAllScores(dataType,grsType,population,'NoCorr')
    tableEURCorr=tableFinalScores=saveAllScores(dataType,grsType,population,'Corr')





