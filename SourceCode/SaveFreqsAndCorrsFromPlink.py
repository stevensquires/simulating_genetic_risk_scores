import pandas_plink
import pandas as pd
import numpy as np
from scipy.stats import pearsonr

def basePath():
    return ''
def path0():
    return basePath()+'Data/GRSinputFiles/'
def path1():
    return basePath()+'Data/OtherData/'
def path2():
    return basePath()+'Data/LDlinkData/'
def pathOut(grsName,popWanted):
    return basePath()+'Data/InputDataOther/1000G/'+grsName+'/'+popWanted+'/Tables/'
def returnData1(plinkName,popWanted,grsFileName,chrPosOrRSID):
    name=path0()+plinkName
    G=pandas_plink.read_plink1_bin(name+'.bed',name+'.bim',name+'.fam',verbose=False)
    bimFile=pd.read_csv(name+'.bim',header=None,sep='\t').set_index(1)
    vals,snps,samples=G.values,list(G.snp.values),list(G.sample.values)
    popData=pd.read_csv(path1()+'1000GpopulationInfo.tsv',sep='\t',index_col=0)
    grsScoreFile=pd.read_excel(path2()+grsFileName).set_index(chrPosOrRSID)
    listScores=list(grsScoreFile.index)
    indices=[]
    for i,sam1 in enumerate(samples):
        row0=popData.loc[sam1]
        if row0['Superpopulation code']==popWanted:
            indices.append(i) 
    return bimFile,vals,snps,samples,listScores,indices,grsScoreFile
def saveFrequencyTable(vals,indices,grsScoreFile,listScores,snps,bimFile,grsName,popWanted):
    tableFrequencies=pd.DataFrame(columns=['A1','A2','Freq1','Freq2'])
    arrayVals=vals[indices,:]
    for i,snp in enumerate(snps):
        split0=snp.split(':')
        snp0=split0[0]+':'+split0[1]
        if snp0 in listScores:
            row0=grsScoreFile.loc[snp0]
            rsNum=row0['RSID']
        freqB=(np.sum(arrayVals[:,i])/(2*arrayVals.shape[0]))
        freqA=1-freqB
        bimRow=bimFile.loc[snp]
        a1,a2=bimRow[4],bimRow[5]
        tableFrequencies.loc[rsNum]=[a1,a2,freqA,freqB]
    tableFrequencies.to_csv(pathOut(grsName,popWanted)+'tabFreqs'+grsName+popWanted+'.csv')
    return tableFrequencies
def saveCorrelationTable(vals,indices,snps,listScores,grsScoreFile,minCorr,grsName,popWanted):
    tableCorrelations=pd.DataFrame(columns=['Chr','SNP1','SNP2','Correlation'])
    arrayVals=vals[indices,:]
    tempList=[]
    for i,snp in enumerate(snps):
        split0=snp.split(':')
        chr1='chr'+str(split0[0])
        snp0=split0[0]+':'+split0[1]
        if snp0 in listScores:
            row0=grsScoreFile.loc[snp0]
            rsNum=row0['RSID']
        tempList.append([chr1,rsNum])
    iter0=0
    for i in range(len(snps)-1):
        vals0=arrayVals[:,i]
        for j in range(i+1,len(snps)):
            vals1=arrayVals[:,j]
            corr,pval=pearsonr(vals0,vals1)
            if abs(corr)>minCorr:
                temp1=tempList[i]
                temp2=tempList[j]
                if temp1[0]!=temp2[0]:
                    print(temp1,temp2)
                tableCorrelations.loc[iter0]=temp1[0],temp1[1],temp2[1],corr
                iter0+=1
    tableCorrelations.to_csv(pathOut(grsName,popWanted)+'tabCorrelations'+grsName+popWanted+'.csv')
    return tableCorrelations

##### This code uses genotyped data from 1000G and calculates frequencies and correlations for superpopulations.
##### The code can be adapted for other array data

superPops=['AMR','AFR','EAS','EUR','SAS']
grsName='GRS2' ### 
plinkName='grs67' #### name of the plink array data
chrPosOrRSID='ChrPos'
for popWanted in superPops:
    bimFile,vals,snps,samples,listScores,indices,grsScoreFile=returnData1(plinkName,grsName,popWanted,chrPosOrRSID)
    tableFreqs=saveFrequencyTable(vals,indices,grsScoreFile,listScores,snps,bimFile,grsName,popWanted)
    minCorr=0.05
    tableCorrelations=saveCorrelationTable(vals,indices,snps,listScores,grsScoreFile,minCorr,grsName,popWanted)
    


















