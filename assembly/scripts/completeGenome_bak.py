import sys
from Bio import SeqIO
import os

installationDirectory = sys.argv[1]
genomeToComplete = sys.argv[2]
removeIntermediateFiles = sys.argv[3]

for seq_record in SeqIO.parse(genomeToComplete,"fasta"):
    genomeSeq = str(seq_record.seq)

firstPortion = genomeSeq[:20000]
outfile = open("tempFasta.fasta","w")
outfile.write(">firstPortion\n"+firstPortion+"\n")
outfile.close()

logFile = open("completeGenome_logfile.txt","w")


#Chop scaffolds:
scaffolds = SeqIO.to_dict(SeqIO.parse("../2_spadesAssembly/scaffolds.fasta","fasta"))
outfile = open("choppedScaffolds.fasta","w")
scafnum = 0
for scaffold in scaffolds:
    scafSeq = str(scaffolds[scaffold].seq)
    if len(scafSeq)>2500:
        for a in range(0,len(scafSeq)-2000,+2000):
            scafnum+=1
            outfile.write(">Scaf_"+str(scafnum)+"\n"+scafSeq[a:a+2000]+"\n")
        scafnum+=1
        outfile.write(">Scaf_"+str(scafnum)+"\n"+scafSeq[a:]+"\n")

#Chop assembled genome:
scaffolds = SeqIO.to_dict(SeqIO.parse(genomeToComplete,"fasta"))
for scaffold in scaffolds:
    scafSeq = str(scaffolds[scaffold].seq)
    for a in range(0,len(scafSeq)-2000,+2000):
        scafnum+=1
        #outfile.write(">Scaf_"+str(scafnum)+"\n"+scafSeq[a:a+2000]+"\n")
    scafnum+=1
    #outfile.write(">Scaf_"+str(scafnum)+"\n"+scafSeq[a:]+"\n")

outfile.close()
os.system("mv choppedScaffolds.fasta allSequences.fasta")

os.system(installationDirectory+"resources/Ragout/bin/ragout --overwrite completeGenomeRagout_recepie")
if os.path.isfile("./ragout-out/allSequences_scaffolds.fasta")==True:

    for seq_record in SeqIO.parse("./ragout-out/allSequences_scaffolds.fasta","fasta"):
        fivePrimeEnd_scaffold = str(seq_record.seq)
        break

    logFile.write("The five prime end was reconstructed with ragout\nand has a length of "+str(len(fivePrimeEnd_scaffold))+"\n\n")

    os.system("makeblastdb -dbtype nucl -in ./ragout-out/allSequences_scaffolds.fasta")
    os.system("blastn -query sequence_trl.fasta -db ragout-out/allSequences_scaffolds.fasta  -task blastn -outfmt 6 -out outputBlast1a.txt")
    blastFile = open("outputBlast1a.txt")
    line = blastFile.readline().rstrip()
    fields = line.split("\t")


    firstPortionReconstructed = fivePrimeEnd_scaffold[int(fields[8]):]
    logFile.write("The five end has been cut according to the position of the sequence a\n and it is now "+str(len(firstPortionReconstructed))+" bp\n\n")
    blastFile.close()

    outfile = open("tempFirstPortion.fasta","w")
    outfile.write(">tempFirstPortion\n"+firstPortionReconstructed+"\n")
    outfile.close()

    outfile = open("firstGenomePortion.fasta","w")
    outfile.write(">first portion of the assembled genome\n"+genomeSeq[500:700]+"\n")
    outfile.close()

    os.system("makeblastdb -dbtype nucl -in tempFirstPortion.fasta")
    os.system("blastn -query firstGenomePortion.fasta -db tempFirstPortion.fasta -outfmt 6 -task blastn -out outputBlast2a.txt")

    blastFile = open("outputBlast2a.txt")
    line = blastFile.readline().rstrip()
    blastFile.close()
    fields = line.split("\t")
    toGetFromReassemble = int(fields[8])
    newGenome = firstPortionReconstructed[:toGetFromReassemble]+genomeSeq[500:]

    newFasta = open("newGenome1.fasta","w")
    newFasta.write(">newFasta\n"+newGenome+"\n")
    newFasta.close()
    logFile.write("The five end has been joined to the previously\nassembled genome that is now  "+str(len(newGenome))+" bp\n\n")

else:
    logfile.write("Something went wrong with the five end ragout reconstruction\n\n")
    exit()










for seq_record in SeqIO.parse("newGenome1.fasta","fasta"):
    genomeSeq = str(seq_record.seq)

lastPortion = genomeSeq[-20000:]
outfile = open("tempFasta.fasta","w")
outfile.write(">lastPortion\n"+firstPortion+"\n")
outfile.close()

#Chop scaffolds:
scaffolds = SeqIO.to_dict(SeqIO.parse("../2_spadesAssembly/scaffolds.fasta","fasta"))
outfile = open("choppedScaffolds.fasta","w")
scafnum = 0
for scaffold in scaffolds:
    scafSeq = str(scaffolds[scaffold].seq)
    if len(scafSeq)>2500:
        for a in range(0,len(scafSeq)-2000,+2000):
            scafnum+=1
            outfile.write(">Scaf_"+str(scafnum)+"\n"+scafSeq[a:a+2000]+"\n")
        scafnum+=1
        outfile.write(">Scaf_"+str(scafnum)+"\n"+scafSeq[a:]+"\n")

#Chop assembled genome:
scaffolds = SeqIO.to_dict(SeqIO.parse(genomeToComplete,"fasta"))
for scaffold in scaffolds:
    scafSeq = str(scaffolds[scaffold].seq)
    for a in range(0,len(scafSeq)-2000,+2000):
        scafnum+=1
        #outfile.write(">Scaf_"+str(scafnum)+"\n"+scafSeq[a:a+2000]+"\n")
    scafnum+=1
    #outfile.write(">Scaf_"+str(scafnum)+"\n"+scafSeq[a:]+"\n")

outfile.close()
os.system("mv choppedScaffolds.fasta allSequences.fasta")

os.system(installationDirectory+"resources/Ragout/bin/ragout --overwrite completeGenomeRagout_recepie2")
if os.path.isfile("./ragout-out/allSequences_scaffolds.fasta")==True:

    for seq_record in SeqIO.parse("./ragout-out/allSequences_scaffolds.fasta","fasta"):
        threePrimeEnd_scaffold = str(seq_record.seq)
        break

    logFile.write("The three prime end was reconstructed with ragout\nand has a length of "+str(len(threePrimeEnd_scaffold))+"\n\n")


    os.system("makeblastdb -dbtype nucl -in ./ragout-out/allSequences_scaffolds.fasta")
    os.system("blastn -query sequence_irs.fasta -db ragout-out/allSequences_scaffolds.fasta  -task blastn -outfmt 6 -out outputBlast1b.txt")
    blastFile = open("outputBlast1b.txt")
    line = blastFile.readline().rstrip()
    fields = line.split("\t")
    print fields[8],fields[9]


    firstPortionReconstructed = threePrimeEnd_scaffold[:int(fields[8])]


    logFile.write("The five end has been cut according to the position of the sequence a\n and it is now "+str(len(firstPortionReconstructed))+" bp\n\n")

    blastFile.close()

    outfile = open("tempFirstPortion.fasta","w")
    outfile.write(">tempFirstPortion\n"+firstPortionReconstructed+"\n")
    outfile.close()

    outfile = open("firstGenomePortion.fasta","w")
    outfile.write(">first portion of the assembled genome\n"+genomeSeq[-700:-500]+"\n")
    outfile.close()

    os.system("makeblastdb -dbtype nucl -in tempFirstPortion.fasta")
    os.system("blastn -query firstGenomePortion.fasta -db tempFirstPortion.fasta -outfmt 6 -task blastn -out outputBlast2b.txt")

    blastFile = open("outputBlast2b.txt")
    line = blastFile.readline().rstrip()
    blastFile.close()
    fields = line.split("\t")
    toGetFromReassemble = int(fields[9])
    newGenome = genomeSeq[:-500]+firstPortionReconstructed[toGetFromReassemble:]

    newFasta = open("newGenome2.fasta","w")
    newFasta.write(">newFasta\n"+newGenome+"\n")
    newFasta.close()

    logFile.write("The three end has been joined to the previously\nassembled genome that is now  "+str(len(newGenome))+" bp\n\n")

else:
    exit()
