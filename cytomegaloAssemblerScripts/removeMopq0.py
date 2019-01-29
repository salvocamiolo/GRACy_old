import sys
from Bio import SeqIO

fasta = open("unmapped.fasta")
sam = open("join.sam")
fasta_out = open("temp.fasta","w")

#Collect sam reads
for a in range(3):
    sam.readline()

seqInSam = {}
unbigousReads = []
num = 0
while True:
    num += 1
    if num % 2000000 == 0:
        print num
    line = sam.readline()
    if not line:
        break
    field = line.split("\t")
    if not field[0] in seqInSam:
        seqInSam[field[0]] = field[4]
    else:
        seqInSam[field[0]+"_2"]= field[4]


for seq_record in SeqIO.parse("unmapped.fasta","fasta"):
    locus = str(seq_record.id)
    if locus in seqInSam and (locus+"_2") in seqInSam:
        doNothing = True
    else:
        fasta_out.write(">"+locus+"\n"+str(seq_record.seq)+"\n")




