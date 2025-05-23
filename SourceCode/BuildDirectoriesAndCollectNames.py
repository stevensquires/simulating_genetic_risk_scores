import pandas as pd
import os
import numpy as np
def basePath():
    return ''
def pathSaveList():
    return basePath()+'Data/OtherData/'
def pathDirectory(pop,grsName):
    return basePath()+'Data/LDlinkData/'+grsName+'/'+pop+'/'
def returnListSuperpopsAndPops():
    ancestryData=pd.read_csv(pathSaveList()+'1000GpopulationInfo.tsv',sep='\t').set_index('Sample name')
    listSuperPops=list(set(ancestryData['Superpopulation code']))
    listPops=list(set(ancestryData['Population code']))
    popsToRemove=['IBS,MSL',np.nan,'GWF','GWJ','GWW','MKK']### remove pops
    for pop in popsToRemove: 
        listPops.remove(pop) 
    listSuperPops.remove('EUR,AFR') ### remove 1 sample which has two superpops
    listSuperPops.remove(np.nan) 
    return listPops,listSuperPops
def buildDirectories0(wantBuildDirectories0):
    if wantBuildDirectories0:
        #### Builds directories as required
        os.mkdir(basePath()+'Data')
        os.mkdir(basePath()+'Data/OtherData/')
        os.mkdir(basePath()+'Data/LDlinkData/')
        os.mkdir(basePath()+'Data/LDlinkData/GRS2/')
def buildDirectories1(wantBuildDirectories1,grsName):
    if wantBuildDirectories1:    
        #### Returns lists of 1000G populations and super-populations. And saves a csv of the available populations
        listPops,listSuperPops=returnListSuperpopsAndPops() ### Produce list of populations and superpopulations from 1000G data
        table1=pd.DataFrame(columns=['Population'],data=listPops)
        table1.to_csv(pathSaveList()+'listPopulations.csv')
        ### Creates folders for saving the LDlink data
        for name in listPops:
            os.mkdir(pathDirectory(name,grsName))
            os.mkdir(pathDirectory(name,grsName)+'Raw')
            os.mkdir(pathDirectory(name,grsName)+'Tables')
        for name in listSuperPops:
            os.mkdir(pathDirectory(name,grsName))
            os.mkdir(pathDirectory(name,grsName)+'Raw')
            os.mkdir(pathDirectory(name,grsName)+'Tables')
                



### Change basePath() to desired path before running


wantBuildDirectories0=False  ### change to true if need to build first set of directories
buildDirectories0(wantBuildDirectories0)

##### need to download the '1000GpopulationInfo.tsv' from 1000G website into basePath()+'Data/OtherData/' folder
##### before running buildDirectories1 as it uses the '1000GpopulationInfo.tsv' file to extract superpopulations
#### and populations

grsName='GRS2'  ### change as desired to make different directories for different GRSs
wantBuildDirectories1=True  ### change to true if want to build all the population and superpop directories for 1000G
buildDirectories1(wantBuildDirectories1,grsName)






