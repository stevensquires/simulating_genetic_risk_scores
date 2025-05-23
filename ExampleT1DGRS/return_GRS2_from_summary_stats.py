import argparse
import pandas as pd
import GenSimulatedArray
import GenerateGRS2


parser = argparse.ArgumentParser(description='Take frequencies and correlations and return a GRS')
parser.add_argument('-t','--tolerance', help='Tolerance required, default=0.01', required=False)
parser.add_argument('-n','--number', help='Number of samples, default=5000', required=False)
parser.add_argument('-i','--iterLimit', help='Iteration limit, default=10000', required=False)
parser.add_argument('-m','--minFreq', help='Minimum frequency, default=0.005', required=False)
parser.add_argument('-p0','--path0', help='Input path for frequency and correlations, default=here', required=False)
parser.add_argument('-p1','--path1', help='Output path for simulated array and GRS, default=here', required=False)
args = vars(parser.parse_args())
# print(args['tolerance'])

if args['tolerance']==None:
    args['tolerance']=0.03
if args['number']==None:
    args['number']=5000
if args['iterLimit']==None:
    args['iterLimit']=30000
if args['minFreq']==None:
    args['minFreq']=0.0001
if args['path0']==None:
    args['path0']=''
if args['path1']==None:
    args['path1']=''

snpsOutOfHWE=pd.read_csv(args['path0']+'SNPsOutOfHWE.txt',header=None)
deviationSNPs={}
for i,row in snpsOutOfHWE.iterrows():
    deviationSNPs[row[0]]=[row[1],row[2],row[3]]
GenSimulatedArray.genSimArr(args,deviationSNPs)
GenerateGRS2.saveAllScores(args)
