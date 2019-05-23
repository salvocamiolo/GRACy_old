#!/usr/bin/python

import os
import sys
from Bio import SeqIO

reads1 = sys.argv[1]
reads2 = sys.argv[2]


notRefined = open("notRefined.txt","w")


os.system("cp finalScaffold.fasta newReference.fasta")
for seq_record in SeqIO.parse("finalScaffold.fasta","fasta"):
            reference = str(seq_record.seq)


newBits = []
previousEnd = []
previousEnd.append(0)
eof = 0
while True and eof == 0:
    for seq_record in SeqIO.parse("newReference.fasta","fasta"):
            reference = str(seq_record.seq)

    os.system("bwaPE newReference.fasta "+reads1+"_subSample.fastq "+reads2+"_subSample.fastq test 16 0.01")
    os.system("samtools faidx newReference.fasta")
    os.system("samtools mpileup -f newReference.fasta test_sorted.bam >finalPileup.txt")
    os.system("maskLowCoverage.py newReference.fasta finalPileup.txt")
    ranges = open("Nranges_smooth.txt")
    line = ranges.readline().rstrip()
    if not line:
        break
    fields = line.split("\t")

    while int(fields[0]) <= previousEnd[len(previousEnd)-1]+50:
        line = ranges.readline().rstrip()
        if not line:
            eof = 1
            break
        fields = line.split("\t")
        
    previousEnd.append(int(fields[1]))
    
    print "REFINING RANGE ",line
    
    
    print "extractSeqByRange.py newReference.fasta finalScaffold "+str(int(fields[0])-3000)+" "+str(int(fields[0])-1000)+" f"
    print "extractSeqByRange.py newReference.fasta finalScaffold "+str(int(fields[1])+1000)+" "+str(int(fields[1])+3000)+" f"
    
    
    if int(fields[0])>3000 and int(fields[0])< (len(reference) - 3000) and  int(fields[1])>3000 and int(fields[1])< (len(reference) - 3000):
        os.system("extractSeqByRange.py newReference.fasta finalScaffold "+str(int(fields[0])-3000)+" "+str(int(fields[0])-1000)+" f")
        os.system("extractSeqByRange.py newReference.fasta finalScaffold "+str(int(fields[1])+1000)+" "+str(int(fields[1])+3000)+" f")
        os.system(" python joinScaffolds_careful.py join "+reads1+".fastq "+reads2+".fastq finalScaffold_"+str(int(fields[0])-3000)+"_"+str(int(fields[0])-1000)+"_f.txt f finalScaffold_"+str(int(fields[1])+1000)+"_"+str(int(fields[1])+3000)+"_f.txt f >joinScaffold_output.txt")
        

        if os.path.isfile("joined_"+"finalScaffold_"+str(int(fields[0])-3000)+"_"+str(int(fields[0])-1000)+"_f.txt_finalScaffold_"+str(int(fields[1])+1000)+"_"+str(int(fields[1])+3000)+"_f.txt")==True:

            for seq_record in SeqIO.parse("joined_"+"finalScaffold_"+str(int(fields[0])-3000)+"_"+str(int(fields[0])-1000)+"_f.txt_finalScaffold_"+str(int(fields[1])+1000)+"_"+str(int(fields[1])+3000)+"_f.txt","fasta"):
                newBit = str(seq_record.seq)
            if len(newBit) > 2500:
                newReference = reference[:reference.find(newBit[:100])]+newBit+reference[reference.find(newBit[-100:])+100:]
                newRef = open("newReference.fasta","w")
                newRef.write(">finalScaffold\n"+newReference+"\n")
                newRef.close()
                print "Finito ciclo ",line
        else:
            print "RANGE NOT CLOSED ",str(int(fields[0])-3000),str(int(fields[0])-1000)
            notRefined.write(str(int(fields[0])-3000)+"\t"+str(int(fields[0])-1000))


        

        


        
        
        
