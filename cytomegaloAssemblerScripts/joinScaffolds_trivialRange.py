import sys
from Bio import SeqIO
import os

scaffold = sys.argv[1]
start = int(sys.argv[2])
end = int(sys.argv[3])

reads1 = sys.argv[4]
reads2 = sys.argv[5]



sequences = {}
readSeq = []


def fuseSequences2(s1,s2):
    start = open("s1.fasta","w")
    start.write(">s1\n"+s1+"\n")
    start.close()
    toFuse = open("s2.fasta","w")
    toFuse.write(">s2\n"+s2+"\n")
    toFuse.close()
    os.system("makeblastdb -dbtype nucl -in s1.fasta")
    os.system("blastn -query s2.fasta -db s1.fasta -outfmt 6 -task blastn  -dust no -soft_masking false -out outputBlast.txt")
    blastFile = open("outputBlast.txt")

    downstreamAlignment = []
    lastNucl = 0

    while True:
        line = blastFile.readline().rstrip()
        if not line:
            break
        fields = line.split("\t")
        if int(fields[9]) > (len(s1)-35) and ( int(fields[9])-int(fields[8]) )>30 and float(fields[10]) <0.0001:
            downstreamAlignment = fields
            lastNucl = int(fields[9])

    if len(downstreamAlignment) >0:
        newSequence = s1[:int(downstreamAlignment[9])]+s2[int(downstreamAlignment[7]):]
        blastFile.close()
        return newSequence 
    else:
        return ""




for seq_record in SeqIO.parse(scaffold,"fasta"):
    scaffoldSeq = str(seq_record.seq)




#print "Collecting reads 1"
#for seq_record in SeqIO.parse(reads1,"fastq"):
#    if not str(seq_record.id) in readSeq:
#        readSeq[str(seq_record.id)] = str(seq_record.seq)

#print "Collecting reads 2"
#for seq_record in SeqIO.parse(reads1,"fastq"):
#    locus = str(seq_record.id)+"_2"
#    if not locus in readSeq:
#        readSeq[locus] = str(seq_record.seq)




elongedSequence = scaffoldSeq[start-500:start]
terminiSeq = scaffoldSeq[end:end+200]
numCycle = 0
while True:
    bestElongation = 0
    numCycle += 1
    print numCycle
    print elongedSequence[-50:]
    os.system("grep "+elongedSequence[-50:]+" "+reads1+" > extractedReads")
    os.system("grep "+elongedSequence[-50:]+" "+reads2+" >> extractedReads")

    tempFile = open("extractedReads")
    readSeq = []
    while True:
        line = tempFile.readline().rstrip()
        if not line:
            break
        readSeq.append(line)
    tempFile.close()


    for seq in readSeq:
        if elongedSequence[-50:] in seq:
            elongationStart = seq.find(elongedSequence[-50:])
            elongation = len(seq) - elongationStart
            if elongation > bestElongation:
                print "Found sequence",seq
                print "Sequence to elong",elongedSequence[-50:]
                
                toElong = elongedSequence + seq[elongationStart+50:-10]
                bestElongation = elongation


    print toElong
    elongedSequence = toElong

    fusedTermini = fuseSequences2(elongedSequence,terminiSeq)
    print "Fine ciclo"
    


    if not fusedTermini == "":
        js = open("joined_"+str(start)+"_"+str(end),"w")
        js.write(">joined_"+str(start)+"_"+str(end)+"\n"+fusedTermini)
        js.close()
        exit()
    















