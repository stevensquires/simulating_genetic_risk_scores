import pandas as pd
import os

file0=pd.read_excel('T2D_GRS392_SEARCH_hg19.xlsx')
dictAll={}
for i,row in file0.iterrows():
    chr0=row['CHR:POS'].split(':')[0]
    if chr0 in dictAll:
        list0=dictAll[chr0]
        list0.append(row['RSID'])
    else:
        dictAll[chr0]=[row['RSID']]

listDirs=os.listdir()
if 'RequestData' not in listDirs:
    os.mkdir('RequestData')
if 'CorrelationData' not in listDirs:
    os.mkdir('CorrelationData')
if 'FrequencyData' not in listDirs:
    os.mkdir('FrequencyData')
listFreqData=os.listdir('FrequencyData')
listCorrData=os.listdir('CorrelationData')

for chrNum,list0 in dictAll.items():
    with open('RequestData/'+'chr'+str(chrNum)+'.txt', 'w') as f:
        for snpName in list0:
            f.write(snpName+'\n')
    if 'chr'+str(chrNum)+'.txt' not in listFreqData:
        with open('FrequencyData/'+'chr'+str(chrNum)+'.txt', 'w') as f:
            f.write('\n')
    if 'chr'+str(chrNum)+'.txt' not in listCorrData:
        with open('CorrelationData/'+'chr'+str(chrNum)+'.txt', 'w') as f:
            f.write('\n')
    


