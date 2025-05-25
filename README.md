This repository contains the code necessary to calculate genetic risk scores (GRSs) from summary statistics from the paper: "Calculating genetic risk scores directly from summary statistics with an application to type 1 diabetes" by Steven Squires, Michael N. Weedon and Richard A. Oram.

The simplest way to use this code is to run, or adapt, code from the ExampleT1DGRS or ExampleT2DGRS folders.

Examples of running the code for two GRSs for type 1 diabetes and type 2 diabetes are in folders ExampleT1DGRS and ExampleT2DGRS respectively. These start from summary statistics (frequencies and correlations) stored in .csv files. Each folder contains a readme file with notes on how to run the code. The code in the ExampleT1DGRS folder, contains code to run the GRS demonstrated in the paper, runs a GRS which includes interaction terms in addition to the standard linear terms and can be adapted to any similar GRS. The code in the ExampleT2DGRS folder.

For further use of the code all relevant code with notes on its use is in the code "SourceCode" which should be useable for anyone familiar with Python. Code is available for the full pipeline from collection of summary statistics to calculation of the final GRS.

A whole-pipeline example is given in FullPipelineExampleT2DGRS with step-by-step instructions on how to calculate a T2DGRS but requires the request of data from ldlink where correlations and frequencies can be collected.
 SNP arrays (or a real array in csv format).
   - If an alternative to the T1DGRS demonstrated here is desired just need to change the GRS score file. If no interaction terms are required then specify that as an option.

