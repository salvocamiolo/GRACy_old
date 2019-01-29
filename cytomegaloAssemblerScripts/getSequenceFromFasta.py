#!/usr/bin/python
from Bio import SeqIO
import sys

title = sys.argv[2]
outfilename = sys.argv[3]

outfile = open(outfilename,"w")
for seq_record in SeqIO.parse(sys.argv[1],"fasta"):
    if title in str(seq_record.id):
        outfile.write(">"+str(seq_record.id)+"\n"+str(seq_record.seq)+"\n")


outfile.close()
