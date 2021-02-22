#!/usr/bin/python3
import argparse
import os
parser = argparse.ArgumentParser(description="description:merge enhancer/promoter")
parser.add_argument("--cong",help="congifure table,#breed;tissue;enhancer.bed;known-promoter.bed;novel-promoter.bed;H3K4me3.bam;H3K27ac.bam;Input.bam;H3K4me3_depth;H3K27ac_depth;Input_depth")
args=parser.parse_args()
dict1={}
dict2={}
dict3={}
dict4={}
dict5={}
dict6={}
dict7={}
dict8={}
dict9={}
with open(args.cong, 'rt') as file1:
    with open(args.cong, 'rt') as file2:
        with open('H3K4me3', 'wt') as H3K4me3_depth:
            with open('H3K27ac', 'wt') as H3K27ac_depth:
                with open('Input', 'wt') as Input_depth:
                    for line2 in file2:
                        line2=line2.rstrip()
                        list2=line2.split("\t")
                        print(list2[8],file=H3K4me3_depth)
                        print(list2[9],file=H3K27ac_depth)
                        print(list2[10],file=Input_depth)
        os.system("cat H3K4me3 Input > H3K4me3-Input.depth")
        os.system("cat H3K27ac Input > H3K27ac-Input.depth")
        os.system("rm H3K27ac H3K4me3 Input")
with open(args.cong, 'rt') as file1:
    with open(args.cong, 'rt') as file2:
        for line1 in file1:
            line1=line1.rstrip()
            list1=line1.split("\t")
            dict1[list1[1]]= ""
            dict2[list1[1]]= ""
            dict3[list1[1]]= ""
            dict4[list1[1]]= ""
            dict5[list1[1]]= ""
            dict6[list1[1]]= ""
            dict7[list1[1]]= ""
            dict8[list1[1]]= ""
            dict9[list1[1]]= ""
        for line2 in file2:
            line2=line2.rstrip()
            list2=line2.split("\t")
            for g in dict1.keys():
                if g == list2[1]:
                    if dict1[g] == "":
                        dict1[g]=list2[2]
                    else:
                        dict1[g]=dict1[g]+' '+list2[2]
            for g in dict2.keys():
                if g == list2[1]:
                    if dict2[g] == "":
                        dict2[g]=list2[3]
                    else:
                        dict2[g]=dict2[g]+' '+list2[3]
            for g in dict3.keys():
                if g == list2[1]:
                    if dict3[g] == "":
                        dict3[g]=list2[4]
                    else:
                        dict3[g]=dict3[g]+' '+list2[4]
            for g in dict4.keys():
                if g == list2[1]:
                    if dict4[g] == "":
                        dict4[g]=list2[6]
                    else:
                        dict4[g]=dict4[g]+' '+list2[6]
            for g in dict5.keys():
                if g == list2[1]:
                    if dict5[g] == "":
                        dict5[g]=list2[7]
                    else:
                        dict5[g]=dict5[g]+' '+list2[7]
            for g in dict6.keys():
                if g == list2[1]:
                    if dict6[g] == "":
                        dict6[g]='T'+'-'+list2[0]+'_'+list2[1]
                    else:
                        dict6[g]=dict6[g]+' '+'T'+'-'+list2[0]+'_'+list2[1]
            for g in dict7.keys():
                if g == list2[1]:
                    if dict7[g] == "":
                        dict7[g]='C'+'-'+list2[0]+'_'+list2[1]
                    else:
                        dict7[g]=dict7[g]+' '+'C'+'-'+list2[0]+'_'+list2[1]
            for g in dict8.keys():
                if g == list2[1]:
                    if dict8[g] == "":
                        dict8[g]=list2[9]
                    else:
                        dict8[g]=dict8[g]+'\n'+list2[9]
            for g in dict9.keys():
                if g == list2[1]:
                    if dict9[g] == "":
                        dict9[g]=list2[10]
                    else:
                        dict9[g]=dict9[g]+'\n'+list2[10]
        for key in dict1.keys():
            cmd="cat {} > {}-enhancer.cat.bed".format(dict1[key],key)
            os.system(cmd)
            cmd="bedtools sort -i {}-enhancer.cat.bed | bedtools merge -c 4 -o collapse -i - > {}-enhancer.merge.bed".format(key,key)
            os.system(cmd)
            cmd="awk '{{if(($2+$3)%2==0){{print $1,($2+$3)/2-1000,($2+$3)/2+1000,$4}} else if(($2+$3)%2!=0){{print $1,($2+$3)/2-999.5,($2+$3)/2+1000.5,$4}}}}' OFS='\t' {}-enhancer.merge.bed |awk '{{if($2<0){{print $1,0,$3,$4}}else if ($2>=0){{print $1,$2,$3,$4}}}}' OFS='\t' - > {}-enhancer.merge.summit.bed".format(key,key)
            os.system(cmd)
            cmd="multiBamSummary BED-file --BED {}-enhancer.merge.summit.bed --bamfiles {} {} --labels {} {} --outRawCounts {}-enhancer.merge.summit.region.readCounts.tab -out {}-enhancer.merge.summit.readCounts.npz".format(key,dict4[key],dict5[key],dict6[key],dict7[key],key,key)
            os.system(cmd)
            import numpy as np
            with open(key, 'wt') as key_depth:
                print(dict8[key]+'\n'+dict9[key],file=key_depth)
            cmd="awk '!/^#/' {}-enhancer.merge.summit.region.readCounts.tab|awk '{{$1=null;$2=null;$3=null;print $0}}' > {}-enhancer_readCounts.matrix".format(key,key)
            os.system(cmd)
            a=np.loadtxt(key+'-enhancer_readCounts.matrix')
            b=np.loadtxt(key)
            np.savetxt(key+'-enhancer_readcounts-depth.matrix',a/b,'%.4f',delimiter='\t')
            cmd="awk '!/^#/' {}-enhancer.merge.summit.region.readCounts.tab|cut -f1,2,3 - >{}-enhancer_readCounts.cordinate".format(key,key)
            os.system(cmd)
            cmd="head -1 {}-enhancer.merge.summit.region.readCounts.tab > {}-enhancer_readCounts.title".format(key,key)
            os.system(cmd)
            cmd="paste -d '\t' {}-enhancer_readCounts.cordinate {}-enhancer_readcounts-depth.matrix | cat {}-enhancer_readCounts.title - > {}-enhancer_readCounts-depth.txt".format(key,key,key,key)
            os.system(cmd)
            cmd="rm {}-enhancer_readCounts.matrix {}-enhancer_readcounts-depth.matrix {}-enhancer_readCounts.cordinate {}-enhancer_readCounts.title {}".format(key,key,key,key,key)
            os.system(cmd)
        for key in dict2.keys():
            cmd="cat {} > {}-known-promoter.cat.bed".format(dict2[key],key)
            os.system(cmd)
            cmd="bedtools sort -i {}-known-promoter.cat.bed | bedtools merge -c 4 -o collapse -i - > {}-known-promoter.merge.bed".format(key,key)
            os.system(cmd)
            cmd="awk '{{if(($2+$3)%2==0){{print $1,($2+$3)/2-1000,($2+$3)/2+1000,$4}} else if(($2+$3)%2!=0){{print $1,($2+$3)/2-999.5,($2+$3)/2+1000.5,$4}}}}' OFS='\t' {}-known-promoter.merge.bed |awk '{{if($2<0){{print $1,0,$3,$4}}else if ($2>=0){{print $1,$2,$3,$4}}}}' OFS='\t' - > {}-known-promoter.merge.summit.bed".format(key,key)
            os.system(cmd)
        for key in dict3.keys():
            cmd="cat {} > {}-novel-promoter.cat.bed".format(dict3[key],key)
            os.system(cmd)
            cmd="bedtools sort -i {}-novel-promoter.cat.bed | bedtools merge -c 4 -o collapse -i - > {}-novel-promoter.merge.bed".format(key,key)
            os.system(cmd)
            cmd="awk '{{if(($2+$3)%2==0){{print $1,($2+$3)/2-1000,($2+$3)/2+1000,$4}} else if(($2+$3)%2!=0){{print $1,($2+$3)/2-999.5,($2+$3)/2+1000.5,$4}}}}' OFS='\t' {}-novel-promoter.merge.bed |awk '{{if($2<0){{print $1,0,$3,$4}}else if ($2>=0){{print $1,$2,$3,$4}}}}' OFS='\t' - > {}-novel-promoter.merge.summit.bed".format(key,key)
            os.system(cmd)
os.system("cat *-enhancer.merge.bed > enhancer.tissue.cat.bed")
os.system("bedtools sort -i enhancer.tissue.cat.bed | bedtools merge -c 4 -o collapse -i - >enhancer.tissue.merge.bed")
os.system("awk '{if(($2+$3)%2==0){print $1,($2+$3)/2-1000,($2+$3)/2+1000,$4} else if(($2+$3)%2!=0){print $1,($2+$3)/2-999.5,($2+$3)/2+1000.5,$4}}' OFS='\t' enhancer.tissue.merge.bed |awk '{if($2<0){print $1,0,$3,$4}else if ($2>=0){print $1,$2,$3,$4}}' OFS='\t' - > enhancer.tissue.merge.summit.bed")
os.system("cat *-known-promoter.merge.bed > known-promoter.tissue.cat.bed")
os.system("bedtools sort -i known-promoter.tissue.cat.bed | bedtools merge -c 4 -o collapse -i - >known-promoter.tissue.merge.bed")
os.system("awk '{if(($2+$3)%2==0){print $1,($2+$3)/2-1000,($2+$3)/2+1000,$4} else if(($2+$3)%2!=0){print $1,($2+$3)/2-999.5,($2+$3)/2+1000.5,$4}}' OFS='\t' known-promoter.tissue.merge.bed |awk '{if($2<0){print $1,0,$3,$4}else if ($2>=0){print $1,$2,$3,$4}}' OFS='\t' - > known-promoter.tissue.merge.summit.bed")
os.system("cat *-novel-promoter.merge.bed > novel-promoter.tissue.cat.bed")
os.system("bedtools sort -i novel-promoter.tissue.cat.bed | bedtools merge -c 4 -o collapse -i - >novel-promoter.tissue.merge.bed")
os.system("awk '{if(($2+$3)%2==0){print $1,($2+$3)/2-1000,($2+$3)/2+1000,$4} else if(($2+$3)%2!=0){print $1,($2+$3)/2-999.5,($2+$3)/2+1000.5,$4}}' OFS='\t' novel-promoter.tissue.merge.bed |awk '{if($2<0){print $1,0,$3,$4}else if ($2>=0){print $1,$2,$3,$4}}' OFS='\t' - > novel-promoter.tissue.merge.summit.bed")
with open(args.cong, 'rt') as file3:
    IPlabel=""
    INPUTlabel=""
    H3K4me3=""
    H3K27ac=""
    Input=""
    for line in file3:
        line=line.rstrip()
        list=line.split("\t")
        if IPlabel=="":
            IPlabel='T'+'-'+list[0]+'-'+list[1]
        else:
            IPlabel=IPlabel+' '+'T'+'-'+list[0]+'-'+list[1]
        if INPUTlabel=="":
            INPUTlabel='C'+'-'+list[0]+'-'+list[1]
        else:
            INPUTlabel=INPUTlabel+' '+'C'+'-'+list[0]+'-'+list[1]
        if H3K4me3=="":
            H3K4me3=list[5]
        else:
            H3K4me3=H3K4me3+' '+list[5]
        if H3K27ac=="":
            H3K27ac=list[6]
        else:
            H3K27ac=H3K27ac+' '+list[6]
        if Input=="":
            Input=list[7]
        else:
            Input=Input+' '+list[7]
    cmd="multiBamSummary BED-file --BED enhancer.tissue.merge.summit.bed --bamfiles {} {} --labels {} {} --outRawCounts enhancer.tissue.merge.summit.region.readCounts.tab -out enhancer.tissue.merge.summit.readCounts.npz".format(H3K27ac,Input,IPlabel,INPUTlabel)
    os.system(cmd)
    cmd="multiBamSummary BED-file --BED known-promoter.tissue.merge.summit.bed --bamfiles {} {} --labels {} {} --outRawCounts known-promoter.tissue.merge.summit.region.readCounts.tab -out known-promoter.tissue.merge.summit.readCounts.npz".format(H3K4me3,Input,IPlabel,INPUTlabel)
    os.system(cmd)
    cmd="multiBamSummary BED-file --BED novel-promoter.tissue.merge.summit.bed --bamfiles {} {} --labels {} {} --outRawCounts novel-promoter.tissue.merge.summit.region.readCounts.tab -out novel-promoter.tissue.merge.summit.readCounts.npz".format(H3K4me3,Input,IPlabel,INPUTlabel)
    os.system(cmd)
import numpy as np
####enhancer read/depth
os.system("awk '!/^#/' enhancer.tissue.merge.summit.region.readCounts.tab|awk '{$1=null;$2=null;$3=null;print $0}' > enhancer_readCounts.matrix")
a=np.loadtxt('enhancer_readCounts.matrix')
b=np.loadtxt('H3K27ac-Input.depth')
np.savetxt("enhancer_readcounts-depth.matrix",a/b,'%.4f',delimiter='\t')
os.system("awk '!/^#/' enhancer.tissue.merge.summit.region.readCounts.tab|cut -f1,2,3 - >enhancer_readCounts.cordinate")
os.system("head -1 enhancer.tissue.merge.summit.region.readCounts.tab > enhancer_readCounts.title")
os.system("paste -d '\t' enhancer_readCounts.cordinate enhancer_readcounts-depth.matrix | cat enhancer_readCounts.title - > enhancer_readCounts-depth.txt")
os.system("rm enhancer_readCounts.matrix enhancer_readcounts-depth.matrix enhancer_readCounts.cordinate enhancer_readCounts.title")
####known-promoter read/depth
os.system("awk '!/^#/' known-promoter.tissue.merge.summit.region.readCounts.tab|awk '{$1=null;$2=null;$3=null;print $0}' > known-promoter_readCounts.matrix")
a=np.loadtxt('known-promoter_readCounts.matrix')
b=np.loadtxt('H3K4me3-Input.depth')
np.savetxt("known-promoter_readcounts-depth.matrix",a/b,'%.4f',delimiter='\t')
os.system("awk '!/^#/' known-promoter.tissue.merge.summit.region.readCounts.tab|cut -f1,2,3 - >known-promoter_readCounts.cordinate")
os.system("head -1 known-promoter.tissue.merge.summit.region.readCounts.tab > known-promoter_readCounts.title")
os.system("paste -d '\t' known-promoter_readCounts.cordinate known-promoter_readcounts-depth.matrix | cat known-promoter_readCounts.title - > known-promoter_readCounts-depth.txt")
os.system("rm known-promoter_readCounts.matrix known-promoter_readcounts-depth.matrix known-promoter_readCounts.cordinate known-promoter_readCounts.title")
####novel-promoter read/depth
os.system("awk '!/^#/' novel-promoter.tissue.merge.summit.region.readCounts.tab|awk '{$1=null;$2=null;$3=null;print $0}' > novel-promoter_readCounts.matrix")
a=np.loadtxt('novel-promoter_readCounts.matrix')
b=np.loadtxt('H3K4me3-Input.depth')
np.savetxt("novel-promoter_readcounts-depth.matrix",a/b,'%.4f',delimiter='\t')
os.system("awk '!/^#/' novel-promoter.tissue.merge.summit.region.readCounts.tab|cut -f1,2,3 - >novel-promoter_readCounts.cordinate")
os.system("head -1 novel-promoter.tissue.merge.summit.region.readCounts.tab > novel-promoter_readCounts.title")
os.system("paste -d '\t' novel-promoter_readCounts.cordinate novel-promoter_readcounts-depth.matrix | cat novel-promoter_readCounts.title - > novel-promoter_readCounts-depth.txt")
os.system("rm novel-promoter_readCounts.matrix novel-promoter_readcounts-depth.matrix novel-promoter_readCounts.cordinate novel-promoter_readCounts.title")
