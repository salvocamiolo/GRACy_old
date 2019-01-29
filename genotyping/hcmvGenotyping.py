#!/usr/bin/python
installationFolder = "/home3/scc20x/Software/mySoftware/GRACy/genotyping"


from Bio import SeqIO
import sys
import os
import datetime
import random as rd

#Initialize variables
orderedHyperLoci = ["rl5a","rl6","rl12","rl13","ul1","ul9","ul11","ul20","ul73","ul74","ul120","ul139","ul146"]
allGroups = set()

read1 = sys.argv[1]
read2 = sys.argv[2]

numReads = int(sys.argv[3])

dbfile = installationFolder+"/"+sys.argv[4]

NumThreads = sys.argv[5]










#****************** START DEDUPLICATION ***********************
#Preprocessing the reads....
#Removing duplicates with fastUniq
inputPath =  '/'.join(read1.split('/')[0:-1])
fileRoot1 = ((read1.split("/"))[-1])
if fileRoot1[-2:] == "fq":
    fileRoot1 = fileRoot1[:-3]
else:
    fileRoot1 = fileRoot1[:-6]

fileRoot2 = ((read2.split("/"))[-1])
if fileRoot2[-2:] == "fq":
    fileRoot2 = fileRoot2[:-3]
else:
    fileRoot2 = fileRoot2[:-6]


outfile = open(fileRoot1+"_IDCard.txt","w")


fileList = "list"+str(rd.randint(0,1000000))+".txt"

os.system("echo "+read1+" >"+fileList)
os.system("echo "+read2+" >>"+fileList)
print "Performing deduplication...."
dedupFile1 = fileRoot1+"_noDup_1.fq"
dedupFile2 = fileRoot2+"_noDup_2.fq"
os.system("fastuniq -i "+fileList+" -t q -o "+dedupFile1+" -p "+dedupFile2)
print "Done!"


#********************   END  DEDUPLICATION ********************************


#******************* Calculate average coverage for deduplicated reads   **************
print "Performing reference alignment...."
os.system("bowtie2 -1 "+dedupFile1+" -2 "+dedupFile2+" -x "+installationFolder+"/fastaFiles/hcmvReference -S alignmenthsbfy43223.sam >null 2>&1")
print "Done!"
print "Converting sam to bam...."
os.system("samtools view -bS -h alignmenthsbfy43223.sam > alignmenthsbfy43223.bam")
print "Done!"
print "Sorting bam...."
os.system("samtools sort -o alignmenthsbfy43223_sorted.bam alignmenthsbfy43223.bam")
print "Done!"
print "Calculating average coverage...."
os.system("samtools depth  alignmenthsbfy43223_sorted.bam  |  awk '{sum+=$3} END { print sum/NR}' >avCoverage.txt")

avCovFile = open("avCoverage.txt")
avCov = float(avCovFile.readline().rstrip())
avCovFile.close()
detectionTreshold = int(avCov*0.2)  #This may be changed by the user
os.system("rm -f alignmenthsbfy43223.* avCoverage.txt "+fileList)

logFile = open(fileRoot1+"_logFile.txt","w")
now = datetime.datetime.now()
logFile.write("Genotyping sample "+fileRoot1+"  started at "+now.strftime("%H:%M")+"\n")


#Collect kmers from database
geneKmers = {}
kmerdbfile = open(dbfile)
kmerdbfile.readline()
while True:
    line = kmerdbfile.readline().rstrip()
    if not line:
        break
    fields = line.split("\t")
    if not (fields[0],fields[1]) in geneKmers:
        geneKmers[(fields[0],fields[1])] = []
        kmerseqs = fields[2].split(",")
        for item in kmerseqs:
            if not len(item) == 0:
                geneKmers[(fields[0],fields[1])].append(item)

print "Reading fastq file...."
#Get sequences in memory
reads = []
numSeq = 0
overallGeneInfo = {}
print "Collecting sequence from first fastq file...."
for seq_record in SeqIO.parse(read1,"fastq"):
    reads.append(str(seq_record.seq))
    reads.append(str(seq_record.seq.reverse_complement()))
    numSeq +=1
    if numSeq == int(numReads/2):
        break

numSeq = 0
print "Collecting sequence from second fastq file...."
for seq_record in SeqIO.parse(read2,"fastq"):
    reads.append(str(seq_record.seq))
    reads.append(str(seq_record.seq.reverse_complement()))
    numSeq +=1
    if numSeq == int(numReads/2):
        break
    
numSeq = 0
for gene in orderedHyperLoci:
    print "Genotyping gene",gene
    
    #Collect specific kmers for the genotypes of this gene
    specificKmerGroup = {}
    for item in geneKmers:
        if item[0] == gene:
            if not item[1] in specificKmerGroup:
                specificKmerGroup[item[1]] = []
            for seqs in geneKmers[item]:
                specificKmerGroup[item[1]].append(seqs)

    countSeq = {}
    
    totCount = 0
    for gr in specificKmerGroup:
        matchingKmer = set()
        allGroups.add(gr)
        if not gr in countSeq:
            countSeq[gr] = 0
        for totSeq in reads:
            for groupSeq in specificKmerGroup[gr]:
                if groupSeq in totSeq:
                    totCount += 1
                    countSeq[gr] +=1
                    matchingKmer.add(groupSeq)
                    break
        if countSeq[gr]> detectionTreshold:
            countSeq[gr]==0
            totCount = totCount - countSeq[gr]
            if len(matchingKmer)>0:
                print "For group",gr,"there are",len(matchingKmer),"matched kmers over a total of",len(geneKmers[(gene,gr)]),". Number of matches:",countSeq[gr]


    #averageCoverage = {} This was used for the previous statistics (not efficient....)
    averageCoverage2 = {}
    totAverageCoverage = 0.0

    outfile.write(gene)
    for gr in countSeq:
        #if not gr in averageCoverage:
        #    averageCoverage[gr] = float(countSeq[gr])/float(len(specificKmerGroup[gr]))
        #    totAverageCoverage += averageCoverage[gr]
        #    print "Statistics",countSeq[gr],len(specificKmerGroup[gr]),totAverageCoverage
        if not gr in averageCoverage2:
            averageCoverage2[gr] = float(countSeq[gr])/float(totCount)
            if averageCoverage2[gr] >0.02: #Here as treshold may be inserted
                print "Values",float(countSeq[gr]),float(totCount)
                outfile.write("\t"+gr+"\t"+str(averageCoverage2[gr]))

            #print "Statistics2",gr,averageCoverage2[gr],countSeq[gr],totCount

    outfile.write("\n")

    #if not gene in overallGeneInfo:
    #    overallGeneInfo[gene] = {}

    #for gr in countSeq:
    #    if totCount>0:
    #        print "The two statistics",gr, float(averageCoverage[gr])/float(totAverageCoverage),averageCoverage2[gr]
    #    if not gr in overallGeneInfo[gene]:
    #        if totCount>0:
    #            overallGeneInfo[gene][gr] = float(averageCoverage[gr])/float(totAverageCoverage)
    #        else:
    #            overallGeneInfo[gene][gr] = 0







#outfile = open(expName+"_IDCard.txt","w")
#outfile.write("Group\t")
#for a in orderedHyperLoci:
#    outfile.write(a+"\t")
#outfile.write("\n")

#for group in allGroups:
#    outfile.write(group+"\t")
#    for a in orderedHyperLoci:
#        if group in overallGeneInfo[a]:
#            if overallGeneInfo[a][group] >0:
#                outfile.write(str(overallGeneInfo[a][group])+"\t")
#            else:
#                outfile.write("\t") 
#        else:
#            outfile.write("\t")
#    outfile.write("\n")
#outfile.close()

now = datetime.datetime.now()
logFile.write("Genotyping sample "+fileRoot1+"  ended at "+now.strftime("%H:%M")+"\n")
logFile.close()
#outfile = open(expName+"_mostFrequentVar.txt","w")
#outfile.write("Gene\t"+expName+"\n")
#for a in range(2,15,+1):
#    os.system("sed 1d "+expName+"_IDCard.txt |cut -f "+str(a)+ " | sort -r -n | head -1 >mostFrequent ")
#    f = open("mostFrequent")
#    mf = float(f.readline().rstrip())
#    f.close()
#    if len(str(mf)) > 2: #this handles the case for which no match was found for any group
#        outfile.write(orderedHyperLoci[a-2]+"\t"+str(mf)+"\n")
#    else:
#        outfile.write(orderedHyperLoci[a-2]+"\t0\n")

#outfile.close()


outfile.close()








        











