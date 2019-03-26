import sys
import os
from Bio import SeqIO
infile=open("./kmerDB/mainDB_seqs.txt")
outfile = open("mainDB_seqs_filtered.txt","w")

line = infile.readline().rstrip()
outfile.write(line+"\n")


while True:
    line = infile.readline().rstrip()
    if not line:
        break
    fields = line.split("\t")
    geneName = fields[0]
    genotypeName = fields[1]
    sequences = fields[2].split(",")
    tempFasta = open("tempFasta.fasta","w")
    numSeq = 0
    for seq in sequences:
        numSeq+=1
        if len(seq)>2:
            tempFasta.write(">Sequence"+str(numSeq)+"\n"+seq+"\n")
    tempFasta.close()
    os.system("dust tempFasta.fasta 3 >maskedFasta.fasta")

    maskedSequence = SeqIO.to_dict(SeqIO.parse("maskedFasta.fasta","fasta"))
    maskedSequenceList = ""
    numMaskedSequence = 0
    numTotSeq = 0
    for sequence in maskedSequence:
        numTotSeq += 1
        if not "N" in  maskedSequence[sequence].seq:
            maskedSequenceList += str(maskedSequence[sequence].seq)
            maskedSequenceList+=","
        else:
            numMaskedSequence += 1
    
    outfile.write(geneName+"\t"+genotypeName+"\t"+maskedSequenceList+"\n")
    print "For gene",geneName,"and genotype",genotypeName,numMaskedSequence,"hav been masked leaving a total of",numTotSeq-numMaskedSequence
    if (numTotSeq-numMaskedSequence) <2:
        print "\n\nATTENTION",(numTotSeq-numMaskedSequence),"\n\n"
    

outfile.close()
