import sys
sys.path.append('../../Analysis1')
import Produce1000GDistributions1 as Analysis1
import pandas as pd
import os

def basePath():
    return '/slade/home/ss1453/Projects/Other/BackgroundT1DGRS/'
def pathDirectory(pop):
    return basePath()+'Data/LDlinkData/GRS2/'+pop+'/'
def pathSaveList():
    return basePath()+'Data/OtherData/'

listPopulations=list(set(Analysis1.returnTable1().set_index('IID')['Population code'].values))
listPopulations.remove('IBS,MSL')
table1=pd.DataFrame(columns=['Population'],data=listPopulations)
table1.to_csv(pathSaveList()+'listPopulations.csv')

for name in listPopulations:
    os.mkdir(pathDirectory(name))
    os.mkdir(pathDirectory(name)+'Raw')
    os.mkdir(pathDirectory(name)+'Tables')
    





