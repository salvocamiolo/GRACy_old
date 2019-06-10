import os
import sys
import datetime
from Bio import SeqIO

#os.system("mkdir -p scripts")

confFile = open("assembly.conf")
logFile = open("logFile.log","w",buffering=0)


#Reading configuration file
#Read data
projectName = ((confFile.readline().rstrip()).split("\t"))[1]
read1 = ((confFile.readline().rstrip()).split("\t"))[1]
read2 = ((confFile.readline().rstrip()).split("\t"))[1]
read1_toFill = ((confFile.readline().rstrip()).split("\t"))[1]
read2_toFill = ((confFile.readline().rstrip()).split("\t"))[1]
confFile.readline() #Read comment

now = datetime.datetime.now()
logFile.write("Date: "+now.strftime("%Y-%m-%d")+"\n")
logFile.write("Sample name: "+projectName+"\n")
logFile.write("Read1 fastq: "+read1+"\n")
logFile.write("Read2 fastq: "+read2+"\n\n\n")




#***************************************************************
#**************** 1 Reads quality filtering ********************
#***************************************************************


qualityFiltering = ((confFile.readline().rstrip()).split("\t"))[1]
logFile.write("Quality filtering: "+qualityFiltering+"\n")

if qualityFiltering == "yes" or qualityFiltering == "Yes":
    now = datetime.datetime.now()
    logFile.write("Quality filtering started at "+now.strftime("%H:%M")+"\n")

    qualConfFile = open("qualityFiltering.conf","w")
    qualConfFile.write("ProjectName\tallSamples\n")
    qualConfFile.write("Sample_start****************************************"+"\n")
    qualConfFile.write("SampleName\t"+projectName+"\n")
    qualConfFile.write("Read1\t"+read1+"\n")
    qualConfFile.write("Read2\t"+read2+"\n")
    qualConfFile.write("Sample_start***Filtering and trimming options\n")
    for a in range(6):
        qualConfFile.write(confFile.readline())
    qualConfFile.write("SampleEnd*******************************************\n")
    qualConfFile.close()
    os.system("python /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/runQualityFiltering.py qualityFiltering.conf")
    os.system("mkdir -p 1_cleanReads")
    #os.system("mv "+projectName+"_hq_?.fastq ./1_cleanReads/")
    #os.system("rm -rf allSamples")
    #performing normalization
    os.system("mv "+projectName+"_hq_1.fastq ./1_cleanReads/qualityFiltered_1.fq")
    os.system("mv "+projectName+"_hq_2.fastq ./1_cleanReads/qualityFiltered_2.fq")
    

    now = datetime.datetime.now()
    logFile.write("Quality filtering ended at "+now.strftime("%H:%M")+"\n\n")
    os.system("rm -f badReads1.fastq badReads2.fastq *filterStats.txt *singletons.fastq paired_normalized.fq paired.fq qualityFiltering.conf")
 
else:
    for a in range(6):
        confFile.readline()
    if os.path.isdir("./1_cleanReads/") == True:
        if os.path.isfile("./1_cleanReads/"+projectName+"_hq_1.fastq") == False:
            print "File ","./1_cleanReads/",projectName,"_hq_1.fastq does not exist, now exiting...."
            exit()
        if os.path.isfile("./1_cleanReads/"+projectName+"_hq_2.fastq") == False:
            print "File ","./1_cleanReads/",projectName,"_hq_2.fastq does not exist, now exiting...."
            exit()
    else:
        os.system("mkdir 1_cleanReads")
        os.chdir("1_cleanReads")
        os.system("ln -s "+read1+" qualityFiltered_1.fq")
        os.system("ln -s "+read2+" qualityFiltered_2.fq")
        os.chdir("../")







#***************************************************************
#*********************** 2 Denovo assembly *********************
#***************************************************************
confFile.readline() #Read comment

denovoAssembly = ((confFile.readline().rstrip()).split("\t"))[1]
#kmer = ((confFile.readline().rstrip()).split("\t"))[1]
if denovoAssembly == "yes" or denovoAssembly=="Yes":
    now = datetime.datetime.now()
    logFile.write("De novo assembly started at "+now.strftime("%H:%M")+"\n")
    print "\nPerforming denovo assembly......."
    os.system("interleave-reads.py 1_cleanReads/qualityFiltered_1.fq 1_cleanReads/qualityFiltered_2.fq -o paired.fq")
    os.system("normalize-by-median.py  -k 17 -C 200 -M 160e9 -p  -o - paired.fq > paired_normalized.fq")
    os.system("python ~/Software/mySoftware/GRACy/assembly/scripts/splitIntervealed.py paired_normalized.fq")
    os.system("mv newRead_1.fastq ./1_cleanReads/"+projectName+"_hq_1.fastq")
    os.system("mv newRead_2.fastq ./1_cleanReads/"+projectName+"_hq_2.fastq")

    os.system("python /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/getBestAssembly.py ./1_cleanReads/"+projectName+"_hq_1.fastq ./1_cleanReads/"+projectName+"_hq_2.fastq assemblyStatistics.txt")
    os.system("mkdir 2_spadesAssembly")
    os.system("mv scaffolds.fasta ./2_spadesAssembly/")
    #os.system("/home3/scc20x/Software/SPAdes-3.12.0-Linux/bin/spades.py -1 ./1_cleanReads/"+projectName+"_hq_1_subSample.fastq -2 ./1_cleanReads/"+projectName+"_hq_2_subSample.fastq --cov-cutoff auto --careful -k "+kmer+" -o 2_spadesAssembly > /home3/scc20x/null 2>&1")
    os.system("rm numReads.txt subsample_1.fq subsample_2.fq null N50.txt")
    os.system("mv assemblyStatistics.txt ./2_spadesAssembly")


else:
    if os.path.isfile("./2_spadesAssembly/scaffolds.fasta") == False:
        print "You chose not to run assembler but scaffolds.fasta file is not there. Now exiting......"
        logFile.write("ERROR!\n You chose not to run assembler but scaffolds.fasta file is not there\n")
        exit()
if os.path.isfile("./2_spadesAssembly/scaffolds.fasta") == False:
    print "Something went wrong with the assembly. Now exiting......"
    logFile.write("ERROR!\n Something went wrong with the assembly\n")
    exit()
else:
    os.system("rm -rf ./2_spadesAssembly/corrected")
    now = datetime.datetime.now()
    logFile.write("De novo assembly ended at "+now.strftime("%H:%M")+"\n\n")






#***************************************************************
#******************* 3 Scaffold Oriantation ********************
#***************************************************************
confFile.readline() #Read comment
performScaffolding = ((confFile.readline().rstrip()).split("\t"))[1]

if performScaffolding == "yes" or performScaffolding == "Yes":
    now = datetime.datetime.now()
    logFile.write("Scaffolding started at "+now.strftime("%H:%M")+"\n\n")
    os.system("mkdir -p 3_scaffoldsOrientation")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/retrieveNodes.py ./3_scaffoldsOrientation/")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/getSequenceFromFasta.py ./3_scaffoldsOrientation/")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/joinScaffolds_careful.py ./3_scaffoldsOrientation/")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/joinScaffolds_trivial.py ./3_scaffoldsOrientation/")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/joinScaffolds.py ./3_scaffoldsOrientation/")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/biomodule.py ./3_scaffoldsOrientation/")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/sequence_a.fasta ./3_scaffoldsOrientation/")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/sequence_a_prime.fasta ./3_scaffoldsOrientation/")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/createCenterScaffold.py ./3_scaffoldsOrientation/")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/merlinGenome_190000_200000_f.txt ./3_scaffoldsOrientation/")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/extractSeqByRange.py ./3_scaffoldsOrientation/")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/gapPrediction.py ./3_scaffoldsOrientation/")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/hcmv_genomes.fasta* ./3_scaffoldsOrientation/")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/gapfillerlibrary ./3_scaffoldsOrientation/")


    print "\nAligning the contigs to the merlin reference gneome......"
    os.chdir("3_scaffoldsOrientation")
    
    os.system("python createCenterScaffold.py "+projectName)
    os.system("cp newFinalScaffold.fasta longestScaffold.fasta ")
    #os.system("scaffold_builder_v2.py -q ../2_spadesAssembly/scaffolds.fasta -r /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/truncatedReference.fasta -p sb")
    #os.system("head -2 sb_Scaffold.fasta >longestScaffold.fasta")
    os.system("python gapPrediction.py longestScaffold.fasta ../1_cleanReads/qualityFiltered_1.fq ../1_cleanReads/qualityFiltered_2.fq")
    
    finalScaffoldFile = open("finalScaffold.fasta","w")
    finalScaffoldFile.write(">finalScaffold\n")
    gapPresent = 0
    for seq_record in SeqIO.parse("filledGenome.fasta","fasta"):
        finalScaffoldFile.write(str(seq_record.seq)+"\n")
        if "N" in str(seq_record.seq):
            gapPresent = 1
    finalScaffoldFile.close()
    if gapPresent == 1:
        os.system("GapFiller -l gapfillerlibrary -s filledGenome.fasta")
        os.system("mv finalScaffold.fasta finalScaffold.fasta_beforeGapfilling")
        os.system("cp standard_output/standard_output.gapfilled.final.fa ./finalScaffold.fasta")
    
    now = datetime.datetime.now()
    logFile.write("Scaffolding finishes at "+now.strftime("%H:%M")+"\n")
    os.chdir("../")
else:
    if os.path.isfile("./3_scaffoldsOrientation/finalScaffold.fasta") == False:
        print "You chose not to run the scaffolding step but finalScaffold.fasta file is not there. Now exiting......"
        logFile.write("You chose not to run the scaffolding step but finalScaffold.fasta file is not there. Now exiting......")
        exit()



#***************************************************************
#********************** 4 Create Consensus *********************
#***************************************************************
confFile.readline() #Read comment
perform1stConsensusCalling = ((confFile.readline().rstrip()).split("\t"))[1]
if perform1stConsensusCalling == "yes" or perform1stConsensusCalling == "Yes":
    now = datetime.datetime.now()
    logFile.write("First consensus calling started at "+now.strftime("%H:%M")+"\n\n")
    os.system("mkdir 4_createConsensus")
    os.system("cp ./3_scaffoldsOrientation/finalScaffold.fasta ./4_createConsensus")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/joinConsensus.py ./4_createConsensus")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/hcmvConsensusCallPipeline ./4_createConsensus")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/bowtiePE ./4_createConsensus")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/createConsensus ./4_createConsensus")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/extractSeqByRange.py ./4_createConsensus")
    os.chdir("4_createConsensus")
    os.system("./hcmvConsensusCallPipeline ../1_cleanReads/qualityFiltered_1.fq ../1_cleanReads/qualityFiltered_2.fq")
    os.system("python joinConsensus.py "+projectName)
    #os.system("./bowtiePE ../3_scaffoldsOrientation/finalScaffold.fasta ../1_cleanReads/"+projectName+"_hq_1.fastq ../1_cleanReads/"+projectName+"_hq_2.fastq test")
    #os.system("./createConsensus ../3_scaffoldsOrientation/finalScaffold.fasta test_sorted.bam")
    #os.system("mv ../3_scaffoldsOrientation/finalScaffold.fasta_con.fasta ./"+projectName+"_genome.fasta")
    now = datetime.datetime.now()
    logFile.write("First consensus calline ended at "+now.strftime("%H:%M")+"\n\n")
    os.chdir("../")
else:
    if os.path.isfile("./4_createConsensus/"+projectName+"_genome.fasta")==False:
        print "You chose not to run the consensus call step but "+projectName+"_genome.fasta  file is not there. Now exiting......"
        logFile.write("You chose not to run the consensus call step but "+projectName+"_genome.fasta  file is not there. Now exiting......")
        exit()


#***************************************************************
#********************** 5 Refine Consensus *********************
#***************************************************************
confFile.readline() #Read comment
refineAssembly = ((confFile.readline().rstrip()).split("\t"))[1]
if refineAssembly == "yes" or refineAssembly == "Yes":
    now = datetime.datetime.now()
    logFile.write("Assembly sequence refining started at "+now.strftime("%H:%M")+"\n\n")

    os.system("mkdir -p 5_refineAssembly")
    os.system("cp ./4_createConsensus/"+projectName+"_genome.fasta ./5_refineAssembly/finalScaffold.fasta")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/joinScaffolds_careful.py ./5_refineAssembly/")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/joinScaffolds.py ./5_refineAssembly/")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/joinScaffolds_trivial.py ./5_refineAssembly/")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/refineAssembly.py ./5_refineAssembly/")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/biomodule.py ./5_refineAssembly/")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/maskLowCoverage.py ./5_refineAssembly/")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/bowtiePE ./5_refineAssembly/")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/bwaPE ./5_refineAssembly/")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/cleanSoftAndUnmapped.py ./5_refineAssembly/")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/sequence_a.fasta ./5_refineAssembly/")
    os.chdir("5_refineAssembly")
    os.system("python refineAssembly.py ../1_cleanReads/qualityFiltered_1.fq ../1_cleanReads/qualityFiltered_2.fq")

    noRef = []
    infile = open("notRefined.txt")
    while True:
        line = infile.readline()
        if not line:
            break
        noRef.append(line)
    infile.close()
    if len(noRef) >0:
        logFile.write("WARNING! The following ranges were not closed during the assembly refining step:\n")
        for line in noRef:
            logFile.write(line)

    for seq_record in SeqIO.parse("newReference.fasta","fasta"):
        newRefSeq = str(seq_record.seq)
    finalScaffoldFile = open("finalScaffold.fasta","w")
    finalScaffoldFile.write(">finalScaffold\n"+newRefSeq+"\n")
    finalScaffoldFile.close()
    

    now = datetime.datetime.now()
    logFile.write("Assembly sequence refining ended at "+now.strftime("%H:%M")+"\n")

    os.chdir("../")

else:
    if os.path.isfile("./5_refineAssembly/finalScaffold.fasta")==False:
        print "You chose not to run the assembly refining step but the finalScaffold.fasta file is not there. Now exiting......"
        logFile.write("You chose not to run the assembly refining step but the finalScaffold.fasta file is not there. Now exiting......")
        exit()


#***************************************************************
#********************** 6 Create  2nd Consensus ****************
#***************************************************************
confFile.readline() #Read comment
perform2ndConsensusCalling = ((confFile.readline().rstrip()).split("\t"))[1]
if perform2ndConsensusCalling == "yes" or perform2ndConsensusCalling == "Yes":
    now = datetime.datetime.now()
    logFile.write("Second consensus calling started at "+now.strftime("%H:%M")+"\n\n")
    os.system("mkdir 6_createConsensus")
    os.system("cp ./5_refineAssembly/finalScaffold.fasta ./6_createConsensus")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/joinConsensus.py ./6_createConsensus")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/hcmvConsensusCallPipeline ./6_createConsensus")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/bowtiePE ./6_createConsensus")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/createConsensus ./6_createConsensus")
    os.system("cp /home3/scc20x/Software/mySoftware/GRACy/assembly/scripts/extractSeqByRange.py ./6_createConsensus")
    os.chdir("6_createConsensus")
    os.system("./hcmvConsensusCallPipeline ../1_cleanReads/qualityFiltered_1.fq ../1_cleanReads/qualityFiltered_2.fq")
    os.system("python joinConsensus.py "+projectName)
    #os.system("./bowtiePE ../5_refineAssembly/finalScaffold.fasta ../1_cleanReads/"+projectName+"_hq_1.fastq ../1_cleanReads/"+projectName+"_hq_2.fastq test")
    #os.system("./createConsensus ../5_refineAssembly/finalScaffold.fasta test_sorted.bam")
    #os.system("mv ../5_refineAssembly/finalScaffold.fasta_con.fasta ./"+projectName+"_genome.fasta")
    now = datetime.datetime.now()
    logFile.write("Second consensus calline ended at "+now.strftime("%H:%M")+"\n\n")
    os.chdir("../")


print "\nPipeline finished!\n"
now = datetime.datetime.now()
logFile.write("Pipeline finishes at "+now.strftime("%H:%M")+"\n")












