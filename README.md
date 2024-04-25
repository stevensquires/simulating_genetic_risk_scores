The code provided here enables the entire pipeline from collection of frequencies and correlations from LDlink or via the genomic data to the generation of the GRS. 

For the 1000G data the code requires a txt file specifying the sample names and populations which can be downloaded from https://www.internationalgenome.org/data-portal/sample (true as of 25/04/2024, website and links may change with time).

For all simulation code of the T1DGRS it requires the T1DGRS score file. There are *** necessary columns: COMPONENT, RSID, SCORE AND SCORE_ALLELE. If format being used is different from RSID (say chr:pos) then that can replace the RSID. To use a different GRS then it can just be directly replaced and should work appropriately. If there are no interaction terms then the COMPONENT column can be excluded.

The code is provided as a set of modules (.py scripts) which can be used individually but also call to one another. There are three main steps:
1) Collection of input data.
2) Generation of simulated SNP arrays.
3) Generation of simulated GRS.

1) Collection of input data.
We provide two ways of collecting the input data required:
a) Via application to LDlink (link works as of 15/04/2024 https://ldlink.nih.gov/?tab=home)
b) Via use of genotype array, which requires the python package pandas_plink (https://pandas-plink.readthedocs.io/en/latest/) to be installed.
1)a) Application to LDlink:
   - If necessary a module allows for automated building of directories and collection of all 26 1000G population names: "BuildDirectoriesAndCollectNames.py"
   - two python modules are "BuildSubmitFrequnciesRequest.py" and "BuildSubmitCorrelationRequest.py" which generate a set of requests for which populations or super-populations are desired and which SNPs are wanted. 
   - two python modules are can then be run to convert the returned data from ldlink into csv tables: "SaveFreqsAsTable.py" and "SaveCorrelationsAsTable.py"

1)b) Use of genotyped array
   - This assumes that the genotyped array is stored as a PLINK bed/bim/fam format.
   - This takes correlations (as Pearson correlations) that are larger than some chosen value - needs to be selected.
   - run "SaveFreqsAndCorrsFromPlink.py" and it will save csv files in the same format as the ldlink produced csv tables.

2) Generation of simulated SNP arrays
   - The frequencies and correlations are saved in csv format. If there are any SNPs which deviate from HWE then they need to be specified in a dict (examples given in code). If no
   - Either "GenSimArrayCorrelations.py" or "GenSimArrayNoCorrelations.py" should be run. Both call to "SimulationAlgorithm.py". 

3) Generation of T1DGRS
   - Run "GenerateGRS2.py" with whatever changes necessary to specify the correct simulated SNP arrays (or a real array in csv format).
   - If an alternative to the T1DGRS demonstrated here is desired just need to change the GRS score file. If no interaction terms are required then specify that as an option.

