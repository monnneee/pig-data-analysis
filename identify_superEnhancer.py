#!/usr/bin/python3
import argparse
import os
parser = argparse.ArgumentParser(description="description:identify superenhancer")
parser.add_argument("--cong",help="congifure table") #configure file is the same as CRE_merge.py
parser.add_argument("--genome",help="SS10")
args=parser.parse_args()
dict1={}
dict2={}
with open(args.cong, 'rt') as file1:
    for line1 in file1:
        line1=line1.rstrip()
        list1=line1.split("\t")
        dict1[list1[0]+'-'+list1[1]]= list1[2]
        dict2[list1[0]+'-'+list1[1]]= list1[6]+"\t"+list1[7]
    for key in dict1.keys():
        cmd="awk '{{print $1"'"\t"'"NR"'"\t"'""'"enhancer"'""'"\t"'"$2"'"\t"'"$3"'"\t"'"$4"'"\t"'""'"."'""'"\t"'""'"."'""'"\t"'"NR}}' {} > {}-enhancer.gff".format(dict1[key],key)
        os.system(cmd)
    for key,value in dict2.items():
        value=value.rstrip()
        list2=value.split("\t")
        cmd="python2.7 /rose/ROSE_main.py -g {} -i {}-enhancer.gff -r {} -o superenhancer -c {}".format(args.genome,key,list2[0],list2[1])
        os.system(cmd)
