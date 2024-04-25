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
def returnFreqNames(grsName,popWanted):
    list0=os.listdir(path0(grsName,popWanted))
    listChrs=[]
    for name in list0:
        if fnmatch(name,'*Frequencies*'):
            listChrs.append(name.split('Frequencies')[0])
    return list(set(listChrs))
def saveFreqTable(grsName,popWanted):
    tableFrequencies=pd.DataFrame(columns=['A1','A2','Freq1','Freq2'])
    listNames=returnFreqNames(grsName,popWanted)
    for chrName in listNames:
        file0=pd.read_csv(path0(grsName,popWanted)+chrName+'Frequencies.txt',sep='\t')
        for i,row in file0.iterrows():
            if row['RS_Number'][0:2]!='rs':
                break
            alleleFreqs=row['Allele Frequency'].split(', ')
            a1,freq1=alleleFreqs[0].split('=')
            a2,freq2=alleleFreqs[1].split('=')
            tableFrequencies.loc[row['RS_Number']]=[a1,a2,freq1,freq2]
    tableFrequencies.to_csv(pathOut(grsName,popWanted)+'tabFreqs'+grsName+popWanted+'.csv')
    


popsWanted=['AFR','SAS','AMR','EAS']
for popWanted in popsWanted:
    grsName='GRS2'
    saveFreqTable(grsName,popWanted)

listPops=pd.read_csv(path1()+'listPopulations.csv')['Population'].tolist()
for popWanted in listPops:
    grsName='GRS2'
    saveFreqTable(grsName,popWanted)






