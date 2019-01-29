#takes in input two fastq files that has been filtered using NGS QC Toolkit each as single end
#and create two fastq files with the mates written in the same order and a third fastq file with the 
#sigletons


import sys
import os
from Bio import SeqIO

filename1 = sys.argv[1]
filename2 = sys.argv[2]

outfileName1 = filename1+"_hq"
outfileName2 = filename2+"_hq"
outfileName3 = filename1+"_singletons.fastq"





loci1 = {}
loci2 = {}

loci1Set = set()
loci2Set = set()


#Load in memory sequences from first fastq file
print "Reading file ",filename1," for sorting"
for seq_record in SeqIO.parse(filename1,"fastq"):
    locusRoot = ((str(seq_record.id)).split("/"))[0]
    if not locusRoot in loci1:
        loci1[locusRoot] = seq_record
        loci1Set.add(locusRoot)
#Load in memory sequences from second fastq file
print "Reading file ",filename2," for sorting"
for seq_record in SeqIO.parse(filename2,"fastq"):
    locusRoot = ((str(seq_record.id)).split("/"))[0]
    if not locusRoot in loci2:
        loci2[locusRoot] = seq_record
        loci2Set.add(locusRoot)

#Calculate common and non common loci
lociCommon = loci1Set & loci2Set
peculiar1 = loci1Set - loci2Set
peculiar2 = loci2Set - loci1Set

#Write the sorted first fastq file
print "Writing file ",outfileName1
with open(outfileName1,"w") as handle :
    for locus in lociCommon:
        SeqIO.write(loci1[locus],handle,"fastq")

#Write the sorted second fastq file
print "Writing file ",outfileName2
with open(outfileName2,"w") as handle :
    for locus in lociCommon:
        SeqIO.write(loci2[locus],handle,"fastq")

#write singletons
print "Writing file ",outfileName3
with open(outfileName3,"w") as handle :
    for locus in peculiar1:
        SeqIO.write(loci1[locus],handle,"fastq")
    for locus in peculiar2:
        SeqIO.write(loci2[locus],handle,"fastq")

#Rename temporary files with original file names
comand = " cp " + outfileName1 + " " + filename1
os.system(comand)
comand = " cp " + outfileName2 + " " + filename2
os.system(comand)
    




