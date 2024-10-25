In this example we take a simple approach to get a T2DGRS out. This GRS was originally developed in this paper: Oram, Richard A., et al. "Utility of diabetes type–specific genetic risk scores for the classification of diabetes type among multiethnic youth." Diabetes care 45.5 (2022): 1124-1131.

Here we manually extract the required SNPs. Several of the scripts are altered to make the process as simple as possible.

1) Run 'ExtractSNPs.py' which produces a folder called RequestData/ which contains text files for the snps divided into chromosomes. It also produces empty txt files in two new folders: FrequencyData and CorrelationData (it will not overwrite any).

2) Go through each file in RequestData folder, e.g. RequestData/chr1.txt, RequestData/chr2.txt etc. and copy the SNPs (in rsid format) into (both):
a) https://ldlink.nih.gov/?tab=ldhap then click on link "Download Variant File" and copy into the empty text files in the FrequencyData folder.
b) https://ldlink.nih.gov/?tab=ldmatrix then click on link "Download R2 File" and copy into the empty text files in the CorrelationData folder

3) Make directory "Tables" then run:
a) SaveFreqsAsTable.py
b) SaveCorrelationsAsTable.py
which produce two csv files in folder Tables.

4) Make directory SimData and run GenSimulatedArray.py which generates SNP arrays and stores them as csv files in SimData folder

5) Make directory Outputs and run  GenerateT2DGRS.py which produces outputs in the Outputs folder.
