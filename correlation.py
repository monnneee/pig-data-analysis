#!/usr/bin/python3
import argparse, os, numpy, pandas, re
parser = argparse.ArgumentParser(description="description:calculate the correlation coefficient between two replicates,and merge their *.Align.tag")
parser.add_argument("--binSize",type=int,default=2000,help="Length in bases of the window used to sample the genome")
parser.add_argument("--numberOfProcessors",type=int,default=10,help="Number of processors to use")
parser.add_argument("--minMappingQuality",type=int,default=25,help="Only reads that have a mapping quality score of at least this are considered")
parser.add_argument("--Rep1BamPath",help="Path of replicate1 bam")
parser.add_argument("--Rep2BamPath",help="Path of replicate2 bam")
args=parser.parse_args()
args.Rep1Label=re.split('[.|_]+',(args.Rep1BamPath.split('/')[-1]))[0]
##print (args.Rep1Label)
args.Rep2Label=re.split('[.|_]+',(args.Rep2BamPath.split('/')[-1]))[0]
args.outFileName=args.Rep1Label+'_merge_'+args.Rep2Label
args.Rep1TagPath='/'.join(args.Rep1BamPath.split('/')[:-1])+'/*.nodup.tagAlign.gz'
##print ('/'.join(args.Rep1BamPath.split('/')[:-1])+'/*.nodup.tagAlign.gz')
args.Rep2TagPath='/'.join(args.Rep2BamPath.split('/')[:-1])+'/*.nodup.tagAlign.gz'
cmd="multiBamSummary bins --binSize {} --numberOfProcessors {} --minMappingQuality {} --bamfiles {} {} --labels {} {} --outFileName {}.readCounts.npz --outRawCounts {}.readCounts.tab".format(args.binSize,args.numberOfProcessors,args.minMappingQuality,args.Rep1BamPath,args.Rep2BamPath,args.Rep1Label,args.Rep2Label,args.outFileName,args.outFileName)
os.system(cmd)
count=pandas.read_table(open(args.outFileName+'.readCounts.tab',"rb"),delimiter="\t",skiprows=0,usecols=[3,4])
with open('correlation_summary.txt','a+') as f1:
    print(args.outFileName+'\t'+str(numpy.min(count.corr())[1]),file=f1)
print (numpy.min(count.corr())[1])
