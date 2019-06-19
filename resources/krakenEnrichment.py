from Bio import SeqIO
import sys
import os

filesList = sys.argv[1]
krakenDB_hcmv = sys.argv[2]
krakenDB_allSpecies = sys.argv[3]
numThreads = int(sys.argv[4])
detectionThreshold = float(sys.argv[5])






infile = open(filesList)
while True:
    read1 = infile.readline().rstrip()
    if not read1:
        break
    read2 = infile.readline().rstrip()
    if not read2:
        break
    

    filename1 = (read1.split("/"))[-1]
    filename2 = (read2.split("/"))[-1]


    print "Loading read 1 into memory"
    read_1 = SeqIO.to_dict(SeqIO.parse(read1,"fastq"))
    print "Loading read 2 into memory"
    read_2 = SeqIO.to_dict(SeqIO.parse(read2,"fastq"))

    print "Performing kraken on hcmv database"
    os.system("kraken -db "+krakenDB_hcmv+"  --threads "+str(numThreads)+"  --paired --fastq-input " +read1+" "+read2+">krakenOutput_hcmv")

    print "Retrieve hcmv reads names"
    kf = open("krakenOutput_hcmv")
    hcmvReadsSet =set()
    while True:
        line = kf.readline().rstrip()
        if not line:
            break
        fields = line.split("\t")
        if fields[2]== '10359':
            hcmvReadsSet.add(fields[1])
    
    print "Generating hcmv specific reads fastq files"
    read1Out = open(filename1+"_hcmvSpecificReads_1.fastq","w")
    read2Out = open(filename2+"_hcmvSpecificReads_2.fastq","w")
    for read in hcmvReadsSet:
    	if len(str(  (read_1[read+"/1"]).seq )) == len( (read_1[read+"/1"]).letter_annotations["phred_quality"] ) and len(str(  (read_2[read+"/2"]).seq )) == len( (read_2[read+"/2"]).letter_annotations["phred_quality"] ):
            SeqIO.write(read_1[read+"/1"],read1Out,"fastq")
            SeqIO.write(read_2[read+"/2"],read2Out,"fastq")

    kf.close()
    print "Performing kraken on secondary database"
    os.system("kraken -db "+krakenDB_allSpecies+"  --threads "+str(numThreads)+"  --paired --fastq-input " +read1+" "+read2+">krakenOutput_all")

    print "Retrieve taxonomy statistics"
    kf = open("krakenOutput_all")
    foundTaxa = {}
    while True:
        line = kf.readline().rstrip()
        if not line:
            break
        fields = line.split("\t")
        if not fields[2] == '10359' and not fields[2] == '9606' and not fields[2] == '0':
            if not fields[2] in foundTaxa:
                foundTaxa[fields[2]] = 0
            foundTaxa[fields[2]] +=1
    
    kf.close()
    print "Tresholds ",len(hcmvReadsSet)*detectionThreshold
    outfile = open(filename1+"_krakenEnrichmentResult.txt","w")
    outfile.write("HCMV\t"+str(len(hcmvReadsSet))+"\n")
    for item in foundTaxa:
        if foundTaxa[item]>len(hcmvReadsSet)*detectionThreshold:
            outfile.write(item+"\t"+str(foundTaxa[item])+"\n")
    outfile.close()
    os.system("rm -f krakenOutput_all krakenOutput_hcmv")
    



    



