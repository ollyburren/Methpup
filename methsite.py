#!/usr/bin/python
import sys

#data=open('/stats/xin/foxp3_16oct2013_analysis/sam_bowtie_local_indel/P1_A10.extendedFrags.noHeader.site')
#output=open('/stats/xin/foxp3_16oct2013_analysis/temp/P1_A10.extendedFrags.noHeader.site','w')
#index=open('/home/xin/FOXP3/project6/NGS6_CpG_sites.txt')
#the gene detected in the sam file must contain in CpG_sites.txt

index=open(sys.argv[1])
data=open(sys.argv[2])
output=open(sys.argv[3],'w')

#make a dictionary containing the ref meth site information 
lookup={}
for a in index.readlines() :
    b=a.strip().split("\t")
    lookup.update({b[0]:b[1:]})

for a in data.readlines() :
    b=a.split()
    #output.write(b[0]+"\t")
    c=[int(x) for x in lookup.get(b[0])]
    str=""
    #if the the first and the last CpG sites are in the sequence
    if c[0]-int(b[1])>0 and c[-1]-int(b[1])<=len(b[2]) :
        for cc in c :
        #note that the index here starts at 0:
            str=str+b[2][cc-int(b[1])]
        output.write(b[0]+"\t"+str+"\n")

output.close()
