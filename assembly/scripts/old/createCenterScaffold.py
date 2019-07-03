from Bio import SeqIO
import sys
import os

projectName = sys.argv[1]
installationDirectory = sys.argv[2]


#new bit
os.system(installationDirectory+"resources/bowtie2-build merlinGenome_190000_200000_f.txt reference")
os.system(installationDirectory+"resources/bowtie2 -x reference -1 ../1_cleanReads/"+projectName+"_hq_1.fastq -2 ../1_cleanReads/"+projectName+"_hq_2.fastq -S alignment.sam")
os.system(installationDirectory+"resources/samtools view -h -Sb alignment.sam >alignment.bam")
os.system(installationDirectory+"resources/samtools view -F 12 -b alignment.bam >botMapped.bam")
os.system(installationDirectory+"resources/samtools view -f 8 -F 4 -b alignment.bam >twoMapped.bam")
os.system(installationDirectory+"resources/samtools merge all.bam botMapped.bam oneMapped.bam twoMapped.bam")
os.system(installationDirectory+"resources/bam2fastq -o reads#.fastq all.bam")
os.system(installationDirectory+"resources/SPAdes-3.12.0-Linux/bin/spades.py -1 reads_1.fastq -2 reads_2.fastq  --cov-cutoff auto --careful -k 51,61,71 -o centerScaffold")
os.system(installationDirectory+"resources/scaffold_builder_v2.py -q ./centerScaffold/scaffolds.fasta -r merlinGenome_190000_200000_f.txt -p sb2 ")




#os.system("scaffold_builder_v2.py -q ../2_spadesAssembly/scaffolds.fasta -r merlinGenome_190000_200000_f.txt -p sb2 >null 2>&1")
#Check the sb2_scaffold file for anomalous bases
outfile = open("temp.fasta","w")
for seq_record in SeqIO.parse("sb2_Scaffold.fasta","fasta"):
    newSeq = ""
    sequence = str(seq_record.seq)
    for a in range(len(sequence)):
        if not sequence[a] =="A" and not sequence[a] =="T" and not sequence[a] =="G" and not sequence[a] =="C" and  not sequence[a] =="a" and not sequence[a] =="t" and not sequence[a] =="g" and not sequence[a] =="c":
            newSeq += "N"
        else:
            newSeq += sequence[a]
    outfile.write(">"+str(seq_record.id)+"\n"+newSeq+"\n")
outfile.close()
os.system("mv temp.fasta sb2_Scaffold.fasta")


os.system(installationDirectory+"resources/lastz merlinGenome_190000_200000_f.txt sb2_Scaffold.fasta --format=general:name1,start1,end1,name2,start2,end2,identity,score --ydrop=50000 >lastzOutput.txt")

bestAlignment = ["scaffold",1,1,1,1,0]

infile = open("lastzOutput.txt")
infile.readline()
while True:
    line = infile.readline().rstrip()
    if not line:
        break
    fields = line.split("\t")
    if int(fields[8])>bestAlignment[5]:
        bestAlignment = [fields[3] , int(fields[1]),int(fields[2]), int(fields[4]),int(fields[5]),int(fields[8])]


print bestAlignment
print "The center scaffold has length",str(bestAlignment[2]-bestAlignment[1])
#print "The best alignment is ",bestAlignment
#Add the center repetitive region to the contigs files
print("./extractSeqByRange.py sb2_Scaffold.fasta "+bestAlignment[0]+" "+str(bestAlignment[3])+" "+str(bestAlignment[4])+" f")
os.system(installationDirectory+"resources/extractSeqByRange.py sb2_Scaffold.fasta "+bestAlignment[0]+" "+str(bestAlignment[3])+" "+str(bestAlignment[4])+" f")
#os.system("cat ../2_spadesAssembly/scaffolds.fasta "+bestAlignment[0]+"_"+str(bestAlignment[3])+"_"+str(bestAlignment[4])+"_f.txt >scaffolds.fasta ")

os.system(installationDirectory+"resources/scaffold_builder_v2.py -q ../2_spadesAssembly/scaffolds.fasta -r /home3/scc20x/hcmvReference/hcmv_genome.fasta -p sb >null 2>&1")

sbFile = open("sb_Scaffold.fasta")
sbFile.readline()
mainScaffold = sbFile.readline().rstrip()
print "The main scaffold has length",len(mainScaffold)
sbFile.close()



#Get the central repetitive region into the variable centralScaffold
for seq_record in SeqIO.parse(bestAlignment[0]+"_"+str(bestAlignment[3])+"_"+str(bestAlignment[4])+"_f.txt","fasta"):
    centerScaffold = str(seq_record.seq)

#for seq_record in SeqIO.parse("sb2_Scaffold.fasta","fasta"):
#    print str(seq_record.id)
#    if str(seq_record.id) == "Scaffold_1":
#        centerScaffold = (str(seq_record.seq))[:10000]

print len(centerScaffold)

#Check for the presence of N stretches and, if present, run gapfiller:
if "N" in centerScaffold:
    print "Ns present in the center sequence. Now running gapfiller...."
    gfFile = open("gapfillerlib.txt","w")
    gfFile.write("lib1 bwa ../1_cleanReads/"+projectName+"_hq_1.fastq ../1_cleanReads/"+projectName+"_hq_2.fastq 300 0.25 FR")
    gfFile.close()
    centerScaffoldFile = bestAlignment[0]+"_"+str(bestAlignment[3])+"_"+str(bestAlignment[4])+"_f.txt"
    os.system(installationDirectory+"resources/GapFiller -l gapfillerlib.txt -s "+centerScaffoldFile)
    os.system("cp ./standard_output/standard_output.gapfilled.final.fa "+centerScaffoldFile)

    


#Get the 5' and 3' recombination coordinates for the central repetitive region centralScaffold
recombinationStart = 0
recombinationEnd = 0
foundPoint = False
for a in range(500):
    if not mainScaffold[10000:].find(centerScaffold[a:a+100]) == -1:
        recombinationStart = mainScaffold[10000:].find(centerScaffold[a:a+100])
        start2 = a
        foundPoint = True

for a in range(500):
    if not mainScaffold[10000:].find(centerScaffold[-100-a:-a]) == -1:
        recombinationEnd = mainScaffold[10000:].find(centerScaffold[-100-a:-a])
        end2=a
        foundPoint = True

print "Recombination start and end for central repetitive region are ",recombinationStart+10000,recombinationEnd+10000
if not recombinationStart==0 and not recombinationEnd==0:
    newScaffold = mainScaffold[:recombinationStart+10000]+centerScaffold[start2:len(centerScaffold)-100-end2]+mainScaffold[recombinationEnd+10000:]
else:
    newScaffold = mainScaffold


outfile = open("newFinalScaffold.fasta","w")
outfile.write(">finalScaffold\n"+newScaffold+"\n")
outfile.close()





