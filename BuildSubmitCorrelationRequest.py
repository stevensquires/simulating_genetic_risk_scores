import pandas as pd
import os
def basePath():
    return '/slade/home/ss1453/Projects/Other/BackgroundT1DGRS/'
def path0():
    return basePath()+'Data/LDlinkData/'
def path1():
    return basePath()+'Data/OtherData/'
def pathOut(popWanted):
    return basePath()+'Data/LDlinkData/GRS2/'+popWanted+'/Raw/'
def returnChrDict(grsFile): ### grsFile 
    dict1={}
    for i,row in grsFile.iterrows():
        chrNum=row['ChrPos'].split(':')[0]
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
def buildRequest(list1,chrSNPDict): ### list1 is list of pops or list of super-populations
    for popWanted in list1:
        with open('SubmitRequests/submitCorrelationRequest'+grsName+popWanted+'.sh', 'w') as f:
            f.write('#!/bin/bash\n')
            for chrNum,listSNPs in chrSNPDict.items():
                sentence=start()+middle(listSNPs)+end(popWanted)+pathOut(popWanted)+'chr'+str(chrNum)+'Correlations.txt\'\n'
                f.write(sentence)

os.mkdir('SubmitRequests')
#### the grs score file should be stored in 'Data/LDlinkData/' folder
nameGRS='grsScores.xlsx' 
#### "grsScores.xlsx" should have column names including 'ChrPos' which should look something like: "6:32432332"
### and 'RSID' which should look something like "rs43432154" (made up rsid and chr:pos)
grsFile=pd.read_excel(path0()+nameGRS)
chrSNPDict=returnChrDict(grsFile)
grsName='GRS2' ### this can be changed if using a different GRS
listSuperPops=['AFR','SAS','AMR','EAS','EUR'] ## 
listPops=pd.read_csv(path1()+'listPopulations.csv')['Population'].tolist()

buildRequest(listSuperPops,chrSNPDict)
buildRequest(listPops,chrSNPDict)
#### submission requests are saved in "SubmitRequests" folder in the directory this code is run in












