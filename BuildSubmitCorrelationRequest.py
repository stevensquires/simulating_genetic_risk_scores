import pandas as pd

def basePath():
    return '/slade/home/ss1453/Projects/Other/BackgroundT1DGRS/'
def path0():
    return basePath()+'Data/LDlinkData/'
def path1():
    return basePath()+'Data/OtherData/'
def pathOut(popWanted):
    return basePath()+'Data/LDlinkData/GRS2/'+popWanted+'/Raw/'
def returnChrDict(grsFile):
    dict1={}
    for i,row in grsFile.iterrows():
        chrNum=row['POSITION_DBSNP151'].split(':')[0]
        if chrNum in dict1:
            list1=dict1[chrNum]
            list1.append(row['RSID'])
        else:
            dict1[chrNum]=[row['RSID']]
    return dict1 

def start():
    return 'curl -k -X GET \'https://ldlink.nih.gov/LDlinkRest/ldmatrix?snps='
def middle(listSNPs):
    str1=listSNPs[0]
    if len(listSNPs)>1:
        for name in listSNPs:
            str1=str1+'%0A'+name
    return str1
def end(popWanted):
    return '&pop='+popWanted+'&r2_d=r2&window=500000&genome_build=grch37&token=b143e06c160a\' > \''

nameGRS='T1DGRS67_1000G_hg19_FinalFor1000GSimulation.xlsx'
grsFile=pd.read_excel(path0()+nameGRS)
chrSNPDict=returnChrDict(grsFile)

grsName='GRS2'

popsWanted=['AFR','SAS','AMR','EAS'] ## ,'EUR'

for popWanted in popsWanted:
    with open('SubmitRequests/submitCorrelationRequest'+grsName+popWanted+'.sh', 'w') as f:
        f.write('#!/bin/bash\n')
        for chrNum,listSNPs in chrSNPDict.items():
            sentence=start()+middle(listSNPs)+end(popWanted)+pathOut(popWanted)+'chr'+str(chrNum)+'Correlations.txt\'\n'
            f.write(sentence)

listPops=pd.read_csv(path1()+'listPopulations.csv')['Population'].tolist()

for popWanted in listPops:
    with open('SubmitRequests/submitCorrelationRequest'+grsName+popWanted+'.sh', 'w') as f:
        f.write('#!/bin/bash\n')
        for chrNum,listSNPs in chrSNPDict.items():
            sentence=start()+middle(listSNPs)+end(popWanted)+pathOut(popWanted)+'chr'+str(chrNum)+'Correlations.txt\'\n'
            f.write(sentence)













