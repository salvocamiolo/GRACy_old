#This software elong the sequences in a fasta file by using the alignment
#of reads from the same organism and exploting the unmapped reads whose
#mate is mapped on the sequence.
#Usage:
#python elong.py projectNAme pathToRead1.fastq pathToRead2.fastq sequencesToElong.fasta numCycles


#Version 0.1



import os
import sys
from Bio import SeqIO
import biomodule

# Import pairwise2 module
from Bio import pairwise2
# Import format_alignment method
from Bio.pairwise2 import format_alignment


projectName =sys.argv[1]
read1 = sys.argv[2]
read2 = sys.argv[3]
sequenceToElong = sys.argv[4]
sequenceToElongOrientation = sys.argv[5]
sequenceToReach = sys.argv[6]
sequenceToReachOrientation = sys.argv[7]


def fuseSequences2(s1,s2):
    #print "S1"
    #print "S2"
    #print s1
    #print s2
    #sys.stdin.read(1)
    start = open("s1.fasta","w")
    start.write(">s1\n"+s1+"\n")
    start.close()
    toFuse = open("s2.fasta","w")
    toFuse.write(">s2\n"+s2+"\n")
    toFuse.close()
    os.system("makeblastdb -dbtype nucl -in s1.fasta  >null 2>&1")
    os.system("blastn -query s2.fasta -db s1.fasta -outfmt 6 -dust no -soft_masking false -task blastn -out outputBlast.txt  >null 2>&1")
    blastFile = open("outputBlast.txt")

    downstreamAlignment = []
    lastNucl = 0

    while True:
        line = blastFile.readline().rstrip()
        if not line:
            break
        fields = line.split("\t")
        if int(fields[9]) > lastNucl and float(fields[2]) > 90.0 and ( int(fields[9])-int(fields[8]) )>30 and float(fields[10]) <0.0001: 
            downstreamAlignment = fields
            lastNucl = int(fields[9])

    if len(downstreamAlignment) >0:
        newSequence = s1[:int(downstreamAlignment[9])]+s2[int(downstreamAlignment[7]):]
        blastFile.close()
        return newSequence 
    else:
	blastFile.close()
        return ""


def greedyElongation(seq):
    unmappedSequences = {}
    reference = str(seq)
    #os.system("mkdir tempor")

    os.system("cd-hit-est  -d 0  -i unmapped.fasta -o unmapped_cdhit.fasta  >null 2>&1")

    cdhitfile = open("unmapped_cdhit.fasta.clstr")
    #Select only high representated unmapped reads

    unmapSeq = {}
    for seq_record in SeqIO.parse("unmapped.fasta","fasta"):
        locus = str(seq_record.id)
        if not locus in unmapSeq:
            unmapSeq[locus] = str(seq_record.seq)


    clusters ={}
    numSeqInCluster = {}

    numCluster = 0
    line = cdhitfile.readline().rstrip()
    if not line:
        print "Do nothing"
    else:
        while True:
            clusterName = "cluster"+str(numCluster)
            if not clusterName in clusters:
                clusters[clusterName] = []
                numSeqInCluster[clusterName] = 0
            line = "start"
            while not line[0] == ">":
                line = cdhitfile.readline().rstrip()
                if not line:
                    break
                numSeqInCluster["cluster"+str(numCluster)] += 1
                clusters[clusterName].append(((line.split(">"))[1].split("..."))[0])
                
                if line[0] == ">":
                    numCluster +=1
            if not line:
                break


    biggerCluster = ""
    seqInBiggerScaffold = 0
    for item in numSeqInCluster:
        if numSeqInCluster[item] > seqInBiggerScaffold:
            seqInBiggerScaffold = numSeqInCluster[item]
            biggerCluster = item
    
    print "Bigger cluster",biggerCluster,"Size",seqInBiggerScaffold
    #sys.stdin.read(1)
    outcdhitfile = open("tempCdhitFile","w")
    if not len(clusters[biggerCluster]) == 0:
        for sequ in clusters[biggerCluster]:
            if sequ in unmapSeq:
                outcdhitfile.write(">"+sequ+"\n"+unmapSeq[sequ]+"\n")
    outcdhitfile.close()
    
    
    #totSeqInClusters = 0
    #for item in numSeqInCluster:
    #    print "Cluster",item,"num",len(clusters[item])
    
    #    if len(clusters[item])>2:
    #        totSeqInClusters += len(clusters[item])
    #        for sequ in clusters[item]:
    #            if sequ in unmapSeq:
    #                outcdhitfile.write(">"+sequ+"\n"+unmapSeq[sequ]+"\n")
    #    if totSeqInClusters >=3:
    #        break
    #outcdhitfile.close()

    os.system("mv tempCdhitFile unmapped.fasta")
    #print "got there"
    #sys.stdin.read(1)





    for seq_record in SeqIO.parse("unmapped.fasta","fasta"):
            locus = str(seq_record.id)
            if not locus in unmappedSequences:
                unmappedSequences[locus] = str(seq_record.seq)
    numElong = 0
    while True:
        #print "Perform the blast of the unmapped sequences"
        toAssemble = open("toAssemble.fasta","w")
        toAssemble.write(">toElong\n"+reference[-200:]+"\n")
        os.system("makeblastdb -in toElong.fasta -dbtype nucl  >null 2>&1")
        os.system("blastn -query unmapped.fasta -db toElong.fasta -outfmt 6 -num_threads 10 -dust no -soft_masking false -out outputBlast.txt >null 2>&1 ")
        #print "done"
        
        #print "Fine Blast"
        #sys.stdin.read(1)
        #print "Fill the toAssemble file"
        blastFile = open("outputBlast.txt")
        while True:
            line = blastFile.readline().rstrip()
            if not line:
                break
            fields = line.split("\t")
            if (int(fields[8]) > (len(reference)-100) or int(fields[9]) > (len(reference)-100)) and abs(int(fields[9])-int(fields[8])) > 40 and float(fields[2])>95:# and fields[5] == "0":

                if int(fields[8]) > int(fields[9]):
                    toAssemble.write(">"+fields[0]+"\n"+biomodule.reverseComplement(unmappedSequences[fields[0]])+"\n")
                else:
                    toAssemble.write(">"+fields[0]+"\n"+unmappedSequences[fields[0]]+"\n")

        toAssemble.close()
        #print "done"
        #os.system("cd-hit-est  -d 0 -c 1 -aS 1 -i toAssemble.fasta -o toAssemble_cdhit.fasta")
        #os.system("mv toAssemble_cdhit.fasta toAssemble.fasta")

        print "Perform phrap assembly step...."
        


        os.system("phrap -minmatch 15  toAssemble.fasta > phrapAssembly 2>null")
        #os.system("cp toAssemble.fasta ./tempor/toAssemble.fasta_"+str(numElong))
        #os.system("cp outputBlast.txt ./tempor/outputBlast.txt_"+str(numElong))
        #os.system("cp toElong.fasta ./tempor/toElong_"+str(numElong))
        numElong += 1

        longestScaffold = ""
        for seq_record in SeqIO.parse("toAssemble.fasta.contigs","fasta"):
            if len(str(seq_record.seq)) >= len(longestScaffold):
                longestScaffold = str(seq_record.seq)

         #Check whether the produced elonged scaffold is on the right orientation
        tempScaffold = open("tempScaffold.fasta","w")
        tempScaffold.write(">tempScaffold\n"+longestScaffold+"\n")
        tempScaffold.close()
        tempScaffold = open("tempScaffold_query.fasta","w")
        tempScaffold.write(">reference\n"+reference[-100:]+"\n")
        tempScaffold.close()
        
        os.system("makeblastdb -in tempScaffold.fasta -dbtype nucl  >null 2>&1")
        os.system("blastn -query tempScaffold_query.fasta -db tempScaffold.fasta -outfmt 6 -num_threads 10 -dust no -soft_masking false -out tempScaffold_outputBlast.txt  >null 2>&1")
        tempScaffold = open("tempScaffold_outputBlast.txt")
        line = tempScaffold.readline().rstrip()
        fieldBlast = line.split("\t")
        if len(fieldBlast)>2:
            if int(fieldBlast[8]) > int(fieldBlast[9]):
                longestScaffold = biomodule.reverseComplement(longestScaffold)
            tempScaffold.close()
            #print "the longest scaffold has size ",len(longestScaffold)
            #print longestScaffold
            #print "done"
            #sys.stdin.read(1)
        else:
            print "WARNING! EXTENSION STOPPED FOR MISSING ELONGMENT!!"
            js = open("joined_W_WARNING_"+sequenceToElong+"_"+sequenceToReach,"w")
            js.write(">joined_W_WARNING_"+sequenceToElong+"_"+sequenceToReach+"\n"+startingSeq[:-1800]+reference)
            js.close()
            exit()


        sc = fuseSequences2(reference,longestScaffold)
        longestScaffold = sc
        print "Elonged sequence has now a size of", len(longestScaffold),"nucleotides" 
        #print longestScaffold
        #sys.stdin.read(1)
        

         #Check whether the produced elonged scaffold is on the right orientation
        #tempScaffold = open("tempScaffold.fasta","w")
        #tempScaffold.write(">tempScaffold\n"+longestScaffold+"\n")
        #tempScaffold.close()
        #os.system("makeblastdb -in tempScaffold.fasta -dbtype nucl")
        #os.system("blastn -query ./startingScaffolds/starting.fasta -db tempScaffold.fasta -dust no -soft_masking false -outfmt 6 -out tempScaffold_outputBlast.txt")
        #tempScaffold = open("tempScaffold_outputBlast.txt")
        #line = tempScaffold.readline().rstrip()
        #if int(fieldBlast[8]) > int(fieldBlast[9]):
        ##fieldBlast = line.split("\t")
        #    longestScaffold = biomodule.reverseComplement(longestScaffold)
        #tempScaffold.close()
        #print "the longest scaffold has size ",len(longestScaffold)
        #print longestScaffold
        #print "done"

        #sys.stdin.read(1)

        
        #print len(longestScaffold),len(reference)
        if len(longestScaffold) <= len(reference):
            #print "WARNING! EXTENSION STOPPED!!"
            #js = open("joined_W_WARNING_"+sequenceToElong+"_"+sequenceToReach,"w")
            #js.write(">joined_W_WARNING_"+sequenceToElong+"_"+sequenceToReach+"\n"+startingSeq[:-1500]+reference)
            #js.close()
            #exit()
            return longestScaffold
        else:
            reference = longestScaffold
            toElong = open("toElong.fasta","w")
            toElong.write(">toElong\n"+reference+"\n")
            toElong.close()
        


sequences = {}
for seq_record in SeqIO.parse(sequenceToElong,"fasta"):
    startingSeq = str(seq_record.seq)
    id1 = str(seq_record.id)
    if sequenceToElongOrientation == "r":
        startingSeq = biomodule.reverseComplement(startingSeq)

for seq_record in SeqIO.parse(sequenceToReach,"fasta"):
    terminiSeq = str(seq_record.seq)
    id2 = str(seq_record.id)
    if sequenceToReachOrientation == "r":
        terminiSeq = biomodule.reverseComplement(terminiSeq)

termfile = open("termini.fasta","w")
termfile.write(">termini\n"+terminiSeq[:500]+"\n")
termfile.close()

toElong = open("toElong.fasta","w")
toElong.write(">toElong\n"+startingSeq[-1800:-300]+"\n")
toElong.close()


numCycle = 0
unmappedContigs = 0
while True:

    for seq_record in SeqIO.parse("toElong.fasta","fasta"):
        reference = str(seq_record.seq)

   # os.system("~/bin/bwa index toElong.fasta")
   # os.system("~/bin/bwa mem toElong.fasta "+ read1 + " " + read2 + " -t 8 > " +projectName+".sam" )
    print "Indexing....."
    os.system("bowtie2-build --threads 10 toElong.fasta toElong  >null 2>&1")
    print "Aligning....."
    os.system(" bowtie2 --local --very-sensitive-local  -p 10 -x toElong -1 "+ read1 + " -2 " + read2 + " -S "+ projectName+".sam >null 2>&1" )
    #os.system("samtools view -bS -h " + projectName+".sam >" + projectName +".bam")
   
    #collect the softclipped reads
    #print "Collecting soft clipped reads"
    alignment = open(projectName+".sam")
    softClipped = open("softclippedReads","w")
    numSoftClipped = 0
    for overlap in range(100,30,-5):
        print "Extending with overalp",overlap
        alignment.seek(0)
        while True:
            line = alignment.readline().rstrip()
            if not line:
                break
            fields = line.split("\t")
            if len(fields)>=6 and not fields[3]=="*":
                if fields[5][-1:] == "S" and fields[5].count('S') == 1 and not "I" in fields[5] and not "D" in fields[5] and not fields[7]=="*":
                    #print "First filter ",fields[9]
                    splitCigar = fields[5].split("M")
                    if len(splitCigar) == 2 and int(splitCigar[0])>overlap  and (int(fields[3])+int(splitCigar[0]))>(len(reference)-3):
                        #print "Second filter ",fields[9]
                        splitCigar2 = splitCigar[1].split("S")
                        if len(splitCigar2) == 2 and int(splitCigar2[0])>25:
                            #print "third filter ",fields[9]
                            #Check the absence of homopolymer at the end of the read
                            if not fields[9][-20:].count("A") > 15 and not fields[9][-20:].count("G") > 15 and not fields[9][-20:].count("C") > 15 and not fields[9][-20:].count("T") > 15:
                                #print "Forth filter ",fields[9]
                                softClipped.write(">"+fields[0]+"_softClipped\n"+fields[9]+"\n")
                                numSoftClipped += 1
                                if numSoftClipped >100:
                                    break
        if numSoftClipped >2:
            break

    softClipped.close()
    os.system("mv  softclippedReads unmapped.fasta")


    longestScaffold = greedyElongation(reference)

    print longestScaffold

    if len(longestScaffold) <= len(reference):
        print "WARNING! EXTENSION STOPPED FOR MISSING ELONGATION!!"
        js = open("joined_W_WARNING_"+sequenceToElong+"_"+sequenceToReach,"w")
        js.write(">joined_W_WARNING_"+sequenceToElong+"_"+sequenceToReach+"\n"+startingSeq[:-1800]+reference)
        js.close()
        exit()

    if len(longestScaffold) > 10000:
        print "WARNING! EXTENSION STOPPED FOR NEVERENDING ELONGATION!!"
        js = open("joined_W_WARNING_"+sequenceToElong+"_"+sequenceToReach,"w")
        js.write(">joined_"+sequenceToElong+"_"+sequenceToReach+"_NeverEndingElongation\n"+startingSeq+"\n")
        js.close()
        exit()

    



    #Check whether the produced elonged scaffold is on the right orientation
    tempScaffold = open("tempScaffold.fasta","w")
    tempScaffold.write(">tempScaffold\n"+longestScaffold+"\n")
    tempScaffold.close()
    os.system("makeblastdb -in tempScaffold.fasta -dbtype nucl >null 2>&1 ")
    os.system("blastn -query toElong.fasta -db tempScaffold.fasta -dust no -soft_masking false -outfmt 6 -out tempScaffold_outputBlast.txt >null 2>&1 ")
    tempScaffold = open("tempScaffold_outputBlast.txt")
    line = tempScaffold.readline().rstrip()
    fieldBlast = line.split("\t")
    if int(fieldBlast[8]) > int(fieldBlast[9]):
        longestScaffold = biomodule.reverseComplement(longestScaffold)
        print "The scaffold is reverse"
    tempScaffold.close()



    #check if the elonged sequence reached the specified termini portion
    #print "Longest scaffold", longestScaffold
    #print "Termini",terminiSeq[:500]
    fusedTermini = fuseSequences2(longestScaffold,terminiSeq[:500])
    #print "Termini checked"
    #sys.stdin.read(1)
    if not fusedTermini == "":
        js = open("joined_"+sequenceToElong+"_"+sequenceToReach,"w")
        js.write(">joined_"+sequenceToElong+"_"+sequenceToReach+"\n"+startingSeq[:-1800]+fuseSequences2(longestScaffold,terminiSeq))
        js.close()
        exit()


    #os.system("cat toElong.fasta termini.fasta > checkTermini.fasta")
    #os.system("phrap -raw checkTermini.fasta")
    #os.system("cp checkTermini.fasta ./tempor/checkTermini.fasta"+str(numCycle))
    #numCycle += 1
    #os.system("mv tempor tempor"+str(numCycle))
    #loci = []
    #for seq_record in SeqIO.parse("checkTermini.fasta.singlets","fasta"):
    #    loci.append(str(seq_record.id))
    #if not "termini" in loci:
    #    endReached = True
    #    js = open("joined_"+id1+"_"+id2,"w")
    #    js.write(">joined_"+id1+"_"+id2+"\n"+fuseSequences2(longestScaffold,terminiSeq))#

    #    js.close()
    #    exit()




 
    





