import pandas as pd
import os
import numpy as np
def basePath():
    return 'C:/Users/ss1453/OneDrive - University of Exeter/Documents/Research/WriteUps/Papers/CurrentlyWorkingOn/SimulatingGRS/TestCodeForGithub/'
def inputPath():
    return basePath()+'InputData/'
def pathSaveList():
    return basePath()+'Data/OtherData/'
def pathDirectory(pop):
    return basePath()+'Data/LDlinkData/GRS2/'+pop+'/'
def returnListSuperpopsAndPops():
    ancestryData=pd.read_csv(inputPath()+'1000GpopulationInfo.tsv',sep='\t').set_index('Sample name')
    listSuperPops=list(set(ancestryData['Superpopulation code']))
    listPops=list(set(ancestryData['Population code']))
    popsToRemove=['IBS,MSL',np.nan,'GWF','GWJ','GWW','MKK']### remove pops
    for pop in popsToRemove: 
        listPops.remove(pop) 
    listSuperPops.remove('EUR,AFR') ### remove 1 sample which has two superpops
    listSuperPops.remove(np.nan) 
    return listPops,listSuperPops
### Change basePath() to desired path before running

#### Builds directories as required
os.mkdir(basePath()+'Data')
os.mkdir(basePath()+'Data/OtherData/')
os.mkdir(basePath()+'Data/LDlinkData/')
os.mkdir(basePath()+'Data/LDlinkData/GRS2/')

#### Returns lists of 1000G populations and super-populations. And saves a csv of the available populations
listPops,listSuperPops=returnListSuperpopsAndPops() ### Produce list of populations and superpopulations from 1000G data
table1=pd.DataFrame(columns=['Population'],data=listPops)
table1.to_csv(pathSaveList()+'listPopulations.csv')
### Creates folders for saving the LDlink data
for name in listPops:
    os.mkdir(pathDirectory(name))
    os.mkdir(pathDirectory(name)+'Raw')
    os.mkdir(pathDirectory(name)+'Tables')
for name in listSuperPops:
    os.mkdir(pathDirectory(name))
    os.mkdir(pathDirectory(name)+'Raw')
    os.mkdir(pathDirectory(name)+'Tables')
        





