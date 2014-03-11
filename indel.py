#!/usr/bin/python

import re
from collections import Counter  #dictionary about how many repeats for each elements---didn't use now
import numpy as np   #cumulative sum
import sys

#data=open('/stats/xin/foxp3_24sep2013_analysis/sam_alignment/P2_F1.extendedFrags.sam')
#output=open('/stats/xin/foxp3_24sep2013_analysis/sam_modify/P2_F1.extendedFrags.sam','w')

#data=open('/stats/xin/foxp3_24sep2013_analysis/sam_bowtie_local/P5_H9.extendedFrags.noHeader.sam')
#output=open('/stats/xin/foxp3_24sep2013_analysis/sam_bowtie_local_indel/P5_H9.extendedFrags.noHeader.sam','w')

data=open(sys.argv[1])
output=open(sys.argv[2],'w')

#a=data.readline()
for a in data.readlines() :
    b=a.split()
    output.write(b[2]+ "\t"+b[3]+"\t")
    #split cigar
    num=re.compile('\d+')
    pattern=re.compile('\D+')
    #numb is the number part in cigar
    numb=[int(y) for y in num.findall(b[5])]
    #patb is the character part in cigar, standing for (mis)match
    patb=pattern.findall(b[5])
    #Check if perfect match
    if len(patb)==1 :
        output.write(b[9]+"\n")
    else :
        seq=list(b[9])
        #qul=list(b[10])
        #check if there is soft clip at the end of the sequence
        if patb[-1]=="S" :
            del seq[-(numb[-1]+1):-1]
        #nd saves the index where there is a deletion in cigar 
        nd=[i for i, j in enumerate(patb) if j == 'D']
        #posb saves the sequence position where there is a (mis)match
        posb=numb[:]
        #check if there is deletion, the number of deletions doesn't count for sequence index
        if len(nd)>0 :
            for i in nd :
                posb[i]=0
        posb=np.cumsum(posb) 
        #check if there is insetion, if there is , replact the insertion with 'I'
        ni=[i for i, j in enumerate(patb) if j == 'I']
        #if there is insertion
        if len(ni)>0 :
            for i in ni :
                seq[posb[i]]='I'
        #check if there is deletion, if there is, make up the deletion with 'D', note that we check from the end to ensure the sequence position information is unchanged before the deletion.
        if len(nd)>0 :
            nd.reverse()
            for i in nd :
                seq.insert(posb[i], 'D'*numb[i])
        #check if there is soft clipping at the beginning
        if patb[0]=="S" :
            seq=seq[posb[0]:]
        #remove when there is an 'I' in seq, and then join seq together as one string
        seq=''.join(list(filter(('I').__ne__, seq)))
        output.write(seq+"\n")

output.close()
