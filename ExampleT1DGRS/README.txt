To run this code with the example frequencies, correlations and SNPs out of HWE run in the command line:

python return_GRS2_from_summary_stats

which produces: SimArr.csv, ScoresPerSNP.csv, FinalScores.npy, FinalScores.csv, DQlistSorted.

Options can be found by running: 

python return_GRS2_from_summary_stats --help

To alter the arguments, for example to change the output path:

python return_GRS2_from_summary_stats -p1 'desired_output_path'

To alter to run this T1DGRS for a different data-set you need to change three files: tabFrequencies.csv; tabCorrelations.csv; SNPsOutOfHWE.txt. 

tabFrequencies.csv: if the SNPs are the same then they do not need to be changed but the frequencies will need changing for the specific data required. It is important to check that A1 (if homozygous on A1 the SNP value is 0) and A2 (if homozygous on A2 the SNP value is 2) are correct as these can be switched even if the SNPs themselves are the same. If the SNPs are different either because of an altered GRS or that alternative SNPs are being used for the same GRS these need to be changed. The SNPs must make with those from the score and interaction files.

tabCorrelations.csv: the pairs of SNPs should be replaced with the correlation between the pair.

SNPsOutOfHWE.txt: this needs to be saved as comma separated txt file with the format: snpName,A1 frequency, heterozygous frequency, A2 frequency.



If you want to alter code to run a linear GRS with no interaction terms use code in ExampleT2DGRS folder.

To alter this code to work for a different GRS with pairwise interaction terms you need to change: 
- grsScores.csv: replace the RSIDs, the SNP weight (SCORE), the effect allele (SCORE_ALLELE) and A1/A2 (with a homozygous A1 having a SNP value of 0 and a homozygous A2 having a SNP value of 2). The interaction term is defined by the COMPONENT column and needs to be linked with interaction_scores.txt and rankingFile.txt.
- An example of an altered GRS is shown with changed grsScores.csv, interaction_scores.txt and rankingFile.txt is shown in folder AlteredFileForAnotherGRS/
