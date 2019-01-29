#!/usr/bin/python

import os
import sys
from Bio import SeqIO
import biomodule

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
    print "Performing Alignment...."
    os.system("./bwaPE newReference.fasta "+reads1+".fastq "+reads2+".fastq test 10 0.02 ")
    os.system("samtools faidx newReference.fasta >null 2>&1")
    print "Calculating coverage on assembly...."
    os.system("samtools mpileup -f newReference.fasta test_sorted.bam >finalPileup.txt 2>null")
    os.system("maskLowCoverage.py newReference.fasta finalPileup.txt ")
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
    
    print "Refining the range",line
    

    
    
    if len(line)>2 and int(fields[0])>3000 and int(fields[0])< (len(reference) - 3000) and  int(fields[1])>3000 and int(fields[1])< (len(reference) - 3000):
        os.system("extractSeqByRange.py newReference.fasta finalScaffold "+str(int(fields[0])-1500)+" "+str(int(fields[0])-500)+" f")
        os.system("extractSeqByRange.py newReference.fasta finalScaffold "+str(int(fields[1])+500)+" "+str(int(fields[1])+1500)+" f")
        
        numTries = 0
        while True:
            numTries += 1
            if numTries ==2:
                break
            print "Performing joinScaffold_careful algorithm on range",line,"...."
            os.system(" python joinScaffolds_careful.py join "+reads1+".fastq "+reads2+".fastq finalScaffold_"+str(int(fields[0])-1500)+"_"+str(int(fields[0])-500)+"_f.txt f finalScaffold_"+str(int(fields[1])+500)+"_"+str(int(fields[1])+1500)+"_f.txt f ")
            if os.path.isfile("joined_finalScaffold_"+str(int(fields[0])-1500)+"_"+str(int(fields[0])-500)+"_f.txt_finalScaffold_"+str(int(fields[1])+500)+"_"+str(int(fields[1])+1500)+"_f.txt")==False:
                print "Refining range",fields[0],fields[1],"with joinScaffold_careful failed. Now trying joinScaffold on range",line,"...."
                os.system("python joinScaffolds.py join "+reads1+".fastq "+reads2+".fastq finalScaffold_"+str(int(fields[0])-1500)+"_"+str(int(fields[0])-500)+"_f.txt f finalScaffold_"+str(int(fields[1])+500)+"_"+str(int(fields[1])+1500)+"_f.txt f")
                if os.path.isfile("joined_finalScaffold_"+str(int(fields[0])-1500)+"_"+str(int(fields[0])-500)+"_f.txt_finalScaffold_"+str(int(fields[1])+500)+"_"+str(int(fields[1])+1500)+"_f.txt")==False:
                    print  "Both alogrithms failed on range,",line
                    print "Now performing 30 steps of joinScaffold trivial in both directions and recording the output"
                    os.system("python joinScaffolds_trivial.py ../1_cleanReads/qualityFiltered_1.fq ../1_cleanReads/qualityFiltered_2.fq finalScaffold_"+str(int(fields[0])-1500)+"_"+str(int(fields[0])-500)+"_f.txt f finalScaffold_"+str(int(fields[1])+500)+"_"+str(int(fields[1])+1500)+"_f.txt f 30" )
                    os.system("mv joinScaffold_trivialSeq.fasta finalScaffold_"+str(int(fields[0])-1500)+"_"+str(int(fields[0])-500)+"_greedy.fasta")
                    os.system("python joinScaffolds_trivial.py ../1_cleanReads/qualityFiltered_1.fq ../1_cleanReads/qualityFiltered_2.fq finalScaffold_"+str(int(fields[1])+500)+"_"+str(int(fields[1])+1500)+"_f.txt r finalScaffold_"+str(int(fields[0])-1500)+"_"+str(int(fields[0])-500)+"_f.txt r 30" )
                    os.system("revcomp joinScaffold_trivialSeq.fasta >temp.fasta; mv temp.fasta joinScaffold_trivialSeq.fasta")
                    os.system("mv joinScaffold_trivialSeq.fasta finalScaffold_"+str(int(fields[1])+500)+"_"+str(int(fields[1])+1500)+"_greedy.fasta")

                else:
                    print "The algorithm joinScaffold was sucessful on range",line
                    break
            else:
                print "The algorithm joinScaffold_careful was sucessful on range",line
                break


        if os.path.isfile("joined_"+"finalScaffold_"+str(int(fields[0])-1500)+"_"+str(int(fields[0])-500)+"_f.txt_finalScaffold_"+str(int(fields[1])+500)+"_"+str(int(fields[1])+1500)+"_f.txt")==True:
            for seq_record in SeqIO.parse("joined_"+"finalScaffold_"+str(int(fields[0])-1500)+"_"+str(int(fields[0])-500)+"_f.txt_finalScaffold_"+str(int(fields[1])+500)+"_"+str(int(fields[1])+1500)+"_f.txt","fasta"):
                newBit = str(seq_record.seq)
            if len(newBit) > 1600:
                newReference = reference[:reference.find(newBit[:100])]+newBit+reference[reference.find(newBit[-100:])+100:]
                newRef = open("newReference.fasta","w")
                newRef.write(">finalScaffold\n"+newReference+"\n")
                newRef.close()

os.system("cp newReference.fasta finalScaffold.fasta")
   

        
        
        
