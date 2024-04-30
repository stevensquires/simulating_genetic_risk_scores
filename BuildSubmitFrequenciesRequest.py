import pandas as pd
def basePath():
    return ''
def path0():
    return basePath()+'Data/LDlinkData/'
def path1():
    return basePath()+'Data/OtherData/'
def path2():
    return basePath()+'Data/OtherData/'
def path3():
    return basePath()+'Code/SimGRS/CollectFreqsAndCorrs/SubmitRequests/'

def pathOut(pop):
    return basePath()+'Data/LDlinkData/GRS2/'+pop+'/Raw/'
def returnChrDict(grsFile):
    dict1={}
    for i,row in grsFile.iterrows():
        chrNum=row['ChrPos'].split(':')[0]
        if chrNum in dict1:
            list1=dict1[chrNum]
            list1.append(row['RSID'])
        else:
            dict1[chrNum]=[row['RSID']]
    #### this below is because LDLink won't accept more than 30 SNPs therefore
    #### split in two/.        
    listChr6=dict1['6']  
    list1=listChr6[:len(listChr6)//2]
    list2=listChr6[len(listChr6)//2:]
    dict1.pop('6')
    dict1['6_1']=list1
    dict1['6_2']=list2
    #####
    return dict1 
def start():
    return 'curl -k -X GET \'https://ldlink.nih.gov/LDlinkRest/ldhap?snps='
def middle(listSNPs):
    str1=listSNPs[0]
    if len(listSNPs)>1:
        for name in listSNPs:
            str1=str1+'%0A'+name
    return str1
def end(popWanted):
    return '&pop='+popWanted+'&r2_d=r2&window=500000&genome_build=grch37&token=b143e06c160a\' > \''
def buildRequest(list1,chrSNPDict):
    for popWanted in popsWanted:
        with open('SubmitRequests/submitFrequencyRequest'+grsName+popWanted+'.sh', 'w') as f:
            f.write('#!/bin/bash\n')
            for chrNum,listSNPs in chrSNPDict.items():
                sentence=start()+middle(listSNPs)+end(popWanted)+pathOut(popWanted)+'chr'+str(chrNum)+'Frequencies.txt\'\n'
                f.write(sentence)



nameGRS='grsScores.xlsx'
grsFile=pd.read_excel(path0()+nameGRS)
chrSNPDict=returnChrDict(grsFile)

grsName='GRS2'
listSuperPops=['AFR','SAS','AMR','EAS','EUR'] ## 
listPops=pd.read_csv(path1()+'listPopulations.csv')['Population'].tolist()
buildRequest(listSuperPops,chrSNPDict)
buildRequest(listPops,chrSNPDict)

#### if below is run then can run this all with one script

with open('SubmitRequests/submitRequestAllPops.sh', 'w') as f:
    f.write('#!/bin/bash\n')
    for popWanted in listPops:
        sentence='bash '+path3()+'submitFrequencyRequestGRS2'+popWanted+'.sh\'\n'
        f.write(sentence)
        sentence='bash '+path3()+'submitCorrelationRequestGRS2'+popWanted+'.sh\'\n'
        f.write(sentence)








