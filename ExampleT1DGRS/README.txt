To run this code with the example frequencies, correlations and SNPs out of HWE run in the command line:

python return_GRS2_from_summary_stats

which produces: SimArr.csv, ScoresPerSNP.csv, FinalScores.npy, FinalScores.csv, DQlistSorted.

Options can be found by running: 

python return_GRS2_from_summary_stats --help

To alter the arguments, for example to change the output path:

python return_GRS2_from_summary_stats -p1 'desired_output_path'

To alter to run this T1DGRS

If you want to alter code to run a linear GRS with no interaction terms use code in ExampleT2DGRS folder.

To alter this code to work for a different GRS with pairwise interaction terms you need to change: 
- grsScores.csv: replace the RSIDs, the SNP weight (SCORE), the effect allele (SCORE_ALLELE) and A1/A2 (with a homozygous A1 having a SNP value of 0 and a homozygous A2 having a SNP value of 2). The interaction term is defined by the COMPONENT column and needs to be linked with interaction_scores.txt and rankingFile.txt
