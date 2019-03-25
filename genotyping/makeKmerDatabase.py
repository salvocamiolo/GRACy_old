from Bio import SeqIO
import sys
import os

orderedHyperLoci = ["rl5a","rl6","rl12","rl13","ul1","ul9","ul11","ul20","ul73","ul74","ul120","ul139","ul146"]

seqFile = sys.argv[1]
kmerSize = int(sys.argv[2])
dbName = sys.argv[3]

numSeq = 0

kmerStatFile = open("./kmerDB/"+dbName+"_stats.txt","w")
kmerStatFile.write("Gene\tGenotype\tNum_of_specific_kmers\n")

kmerSeqFile = open("./kmerDB/"+dbName+"_seqs.txt","w")
kmerSeqFile.write("Gene\tGenotype\tSequences\n")

for gene in orderedHyperLoci:
    print "Creating Kmer database for gene",gene,"...."
    filename = "./fastaFiles/"+gene+"_multiFasta.fasta"
    outfile = open(filename,"w")
    for seq_record in SeqIO.parse(seqFile,"fasta"):
        splitLocusName = str(seq_record.id).split("_")
        if splitLocusName[0] == gene:
            outfile.write(">"+splitLocusName[2]+"_"+str(numSeq)+"\n"+str(seq_record.seq)+"\n")
            numSeq+=1
    outfile.close()
    gKmer = {}
    sequences = {} 
    print "Calculating kmer for all sequences...."
    #Calculate kmer for all sequences
    for seq_record in SeqIO.parse("./fastaFiles/"+gene+"_multiFasta.fasta","fasta"):
        group =(str(seq_record.id).split("_"))[0]

        if not group in gKmer:
            gKmer[group] = set()
        if not group in sequences:
            sequences[group] = []
        seq = str(seq_record.seq)
        sequences[group].append(seq)
        for a in range(0,len(seq)-kmerSize-1,1):
            gKmer[group].add(seq[a:a+kmerSize])
    print "Finding specific kmer per sequence...."
    #Calculate specific kmer for sequence
    specificKmerSequence = {}
    for item in gKmer:
        if not item in specificKmerSequence:
            specificKmerSequence[item] = gKmer[item]
            for item2 in gKmer:
                if not item == item2:
                    specificKmerSequence[item] = specificKmerSequence[item] - gKmer[item2]


    #Calculate specific kmer for group
    print "Finding speicific kmer for group...."
    specificKmerGroup = {}
    for gr in specificKmerSequence:
        if  not gr in specificKmerGroup:
            specificKmerGroup[gr] = set()
        for specificSeq in specificKmerSequence[gr]:
            specific = True
            for sequence in sequences[gr]:
                if not specificSeq in sequence:
                    specific=False
            if specific == True:
                specificKmerGroup[gr].add(specificSeq)

    deletedGroups = []
    for item in specificKmerGroup:
        if not len(specificKmerGroup[item]) > 0:
            print "Group ",item,"has not specific kmer and will be removed from the analysis"
            deletedGroups.append(item)
    for item in deletedGroups:
        del specificKmerGroup[item]



    for item in specificKmerGroup:
        print item,len(specificKmerGroup[item])
        kmerStatFile.write(gene+"\t"+item+"\t"+str(len(specificKmerGroup[item]))+"\n")
        kmerSeqFile.write(gene+"\t"+item+"\t")
        for kmerseq in specificKmerGroup[item]:
            kmerSeqFile.write(kmerseq+",")
        kmerSeqFile.write("\n")

        




