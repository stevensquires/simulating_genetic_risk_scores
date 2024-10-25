import pandas as pd
import os
from fnmatch import fnmatch

def basePath():
    return ''
def path0():
    return basePath()+'FrequencyData/'
def pathOut():
    return basePath()+'Tables/'
def saveFreqTable(grsName,popWanted):
    tableFrequencies=pd.DataFrame(columns=['A1','A2','Freq1','Freq2'])
    listNames=os.listdir(path0())
    for name in listNames:
        file0=pd.read_csv(path0()+name,sep='\t')
        for i,row in file0.iterrows():
            if row['RS_Number'][0:2]!='rs':
                break
            alleleFreqs=row['Allele Frequency'].split(', ')
            a1,freq1=alleleFreqs[0].split('=')
            a2,freq2=alleleFreqs[1].split('=')
            tableFrequencies.loc[row['RS_Number']]=[a1,a2,freq1,freq2]
    tableFrequencies.to_csv(pathOut()+'tabFreqs'+grsName+popWanted+'.csv')
    
grsName='T2DGRS'
popWanted='EUR'
saveFreqTable(grsName,popWanted)







