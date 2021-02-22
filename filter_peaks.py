#!/usr/bin/python3
import argparse
import os,sys
parser = argparse.ArgumentParser(description="description:filter peak by p-value & rpkm,and merge peak")
parser.add_argument("--marker",help="eg.LW-2W-1-H3K4me3-Heart")
parser.add_argument("--PeakPath",help="Path of peak_xls")
parser.add_argument("--IpBamPath",help="Path of Treat bam")
parser.add_argument("--InputbamPath",help="Path of Input bam")
parser.add_argument("--IpDepth",type=float,help="read depth of Treat")
parser.add_argument("--InputDepth",type=float,help="read depth of Input")
parser.add_argument("--distance",type=int,default=0,help="Maximum distance between peaks")
args=parser.parse_args()
## filter peak by p-value & get summit
######
argv0_list = sys.argv[0].split("\\")
script_name = argv0_list[0]
######
cmd="awk '!/^#/ && $7>5' {} |cut -f1,2,3,10| awk '{{if(($2+$3)%2==0){{print $1,($2+$3)/2-1000,($2+$3)/2+1000,$4}} else if(($2+$3)%2!=0){{print $1,($2+$3)/2-999.5,($2+$3)/2+1000.5,$4}}}}' OFS='\t' - |awk '{{if($2<0){{print $1,0,$3,$4}}else if ($2>=0){{print $1,$2,$3,$4}}}}' OFS='\t' - > pFilter-{}.bed".format(args.PeakPath,args.marker)
os.system(cmd)
## reads on peaks
cmd="multiBamSummary BED-file --BED pFilter-{}.bed --bamfiles {} {} --labels Treat Contrl --outRawCounts pFilter-{}.readCounts.tab -out pFilter-{}.readCounts.npz".format(args.marker,args.IpBamPath,args.InputbamPath,args.marker,args.marker)
os.system(cmd)
## match peak with name
cmd="bedtools sort -i pFilter-{}.bed|awk '{{print $1"'":"'"$2"'":"'"$3"'"\t"'"$4}}' - >pFilter-{}.sort.bed".format(args.marker,args.marker)
os.system(cmd)
cmd="bedtools sort -i pFilter-{}.readCounts.tab|awk '!/^#/ {{print $1"'":"'"$2"'":"'"$3"'"\t"'"$4"'"\t"'"$5}}' - >pFilter-{}.readCounts.sort.tab".format(args.marker,args.marker)
os.system(cmd)
cmd="join pFilter-{}.readCounts.sort.tab pFilter-{}.sort.bed |sed 's/:/\t/g' -> pFilter-{}.readCounts.match.tab".format(args.marker,args.marker,args.marker)
os.system(cmd)
## filter peak by rpm
cmd="awk '!/^#/' pFilter-{}.readCounts.match.tab | awk '{{if($5!=0){{print $1,$2,$3,$4,$5,$4/{},$5/{},($4/{})/($5/{}),$4/{}-$5/{},$6}} else if($5==0){{print $1,$2,$3,$4,$5,$4/{},$5/{},"'"NA"'",$4/{}-$5/{},$6}}}}' OFS='\t' - |awk '$8>2 && $9>1' |cut -f1,2,3,10> p-rpkmFilter-{}.bed".format(args.marker,args.IpDepth,args.InputDepth,args.IpDepth,args.InputDepth,args.IpDepth,args.InputDepth,args.IpDepth,args.InputDepth,args.IpDepth,args.InputDepth,args.marker)
os.system(cmd)
##merge peak
cmd="bedtools sort -i p-rpkmFilter-{}.bed | bedtools merge -c 4 -o collapse -d {} -i - > {}.merge.bed".format(args.marker,args.distance,args.marker)
os.system(cmd)
## get summit of merged peak
cmd="awk '{{if(($2+$3)%2==0){{print $1,($2+$3)/2-1000,($2+$3)/2+1000,$4}} else if(($2+$3)%2!=0){{print $1,($2+$3)/2-999.5,($2+$3)/2+1000.5,$4}}}}' OFS='\t' {}.merge.bed| awk '{{if($2<0){{print $1,0,$3,$4}}else if ($2>=0){{print $1,$2,$3,$4}}}}' OFS='\t' - > {}.summit2.bed".format(args.marker,args.marker)
os.system(cmd)
with open(args.marker + ".sh",'wt') as file1:
    print("python3 "+script_name +"  --marker "+ args.marker + " --PeakPath " + args.PeakPath + " --IpBamPath " + args.IpBamPath + " --InputbamPath " + args.InputbamPath + " --IpDepth " + str(args.IpDepth) + " --InputDepth " + str(args.InputDepth),file=file1)
cmd="mkdir {}".format(args.marker)
os.system(cmd)
cmd="rm pFilter-{}.readCounts.npz".format(args.marker)
os.system(cmd)
cmd="mv pFilter-{}.bed pFilter-{}.readCounts.tab pFilter-{}.sort.bed pFilter-{}.readCounts.sort.tab pFilter-{}.readCounts.match.tab p-rpkmFilter-{}.bed {}.merge.bed {}.summit2.bed {}.sh {}".format(args.marker,args.marker,args.marker,args.marker,args.marker,args.marker,args.marker,args.marker,args.marker,args.marker)
os.system(cmd)
OutFile1="{}/p-rpkmFilter-{}.bed".format(args.marker,args.marker)
OutFile2="{}/{}.summit2.bed".format(args.marker,args.marker)
count1=len(open(OutFile1).readlines())
count2=len(open(OutFile2).readlines())
with open('FilterPeak_number.txt','a+') as f1:
    print(os.getcwd()+'\t'+args.marker+'\t'+OutFile1+'\t'+str(count1)+'\t'+OutFile2+'\t'+str(count2),file=f1)
