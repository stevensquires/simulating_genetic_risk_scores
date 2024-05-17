import pandas as pd
import os
def basePath():
    return ''
def path0():
    return basePath()+'Data/LDlinkData/'
def path1():
    return basePath()+'Data/OtherData/'
def pathOut(popWanted,grsName):
    return basePath()+'Data/LDlinkData/'+grsName+'/'+popWanted+'/Raw/'
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
def end(popWanted,token):
    return '&pop='+popWanted+'&r2_d=r2&window=500000&genome_build=grch37&token='+token+'\' > \''
def buildRequest(list1,chrSNPDict,grsName,token): ### list1 is list of pops or list of super-populations
    for popWanted in list1:
        with open(path0()+'SubmitRequests/submitCorrelationRequest'+grsName+popWanted+'.sh', 'w') as f:
            f.write('#!/bin/bash\n')
            for chrNum,listSNPs in chrSNPDict.items():
                sentence=start()+middle(listSNPs)+end(popWanted,token)+pathOut(popWanted,grsName)+'chr'+str(chrNum)+'Correlations.txt\'\n'
                f.write(sentence)
if not os.path.exists(path0()+'SubmitRequests'):
    os.mkdir(path0()+'SubmitRequests')
#### the grs score file should be stored in 'Data/LDlinkData/' folder
nameGRS='grsScores.xlsx' 
#### "grsScores.xlsx" should have column names including 'ChrPos' which should look something like: "6:32432332"
### and 'RSID' which should look something like "rs43432154" (made up rsid and chr:pos) to work for this code
grsFile=pd.read_excel(path0()+nameGRS)
chrSNPDict=returnChrDict(grsFile)
grsName='GRS2' ### this can be changed if using a different GRS
listSuperPops=['AFR','SAS','AMR','EAS','EUR'] ## 
listPops=pd.read_csv(path1()+'listPopulations.csv')['Population'].tolist()

token='b143e06c160a' ### this needs to be generated from the ldlink website to give you access to the data

buildRequest(listSuperPops,chrSNPDict,grsName,token)
buildRequest(listPops,chrSNPDict,grsName,token)
#### submission requests are saved in "SubmitRequests" folder in the directory this code is run in












