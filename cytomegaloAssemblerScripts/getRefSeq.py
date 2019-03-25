import biomodule
from Bio import SeqIO
import sys

contigs = {}

for seq_record in SeqIO.parse("scaffolds.fasta","fasta"):
    locus = str(seq_record.id)
    seq = str(seq_record.seq)
    if not locus in contigs:
        contigs[locus] = seq

infile = open("center.coords")
line = infile.readline()
while not line[:5]=="=====":
    line = infile.readline()

alignments = []

while True:
    line = infile.readline().rstrip()
    newline = line.replace('\t',' ')
    line = newline
    if not line:
        break
    split1 = line.split(" ")
    split2 = []
    for item in split1:
        if not item == '|' and not item == "" and not item == "\t":
            split2.append(item)
    alignments.append(split2)
    
startPos = 0
newSequence = ""

newAlignments =  sorted(alignments, key=lambda tup: tup[0])


for item in newAlignments:
    scarto=0
    if int(item[0])>startPos:
        for a in range(startPos,int(item[0])):
            newSequence+="N"
    else:
        scarto=startPos-int(item[0])
    if int(item[2])<int(item[3]):
        newSequence+=contigs[item[8]][int(item[2])-scarto:int(item[3])]
    else:
        newSequence+=biomodule.reverseComplement(contigs[item[8]][int(item[3]):int(item[2])-scarto])
    startPos = int(item[1])

newSeqFile = open("centerSequence.fasta","w")
newSeqFile.write(">Center_sequence\n"+newSequence+"\n")
newSeqFile.close()

    