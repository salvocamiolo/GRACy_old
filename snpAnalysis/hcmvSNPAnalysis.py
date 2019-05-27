import os
import sys

referenceFile = sys.argv[1] #change in GUI
inputFiles = sys.argv[2]#change in GUI
outputFolder = sys.argv[3] #change in GUI
vcfFileInput = sys.argv[4] #change in GUI
numThreads = sys.argv[5] #change in GUI


if vcfFileInput == '0': #Fastq files are provided and alignments/SNP calling need to be computed
    
    os.system("bowtie2-build "+referenceFile+" reference -q")
    
    infile = open(inputFiles)
    while True:
        read1 = infile.readline().rstrip()
        if not read1:
            break

        read2 = infile.readline().rstrip()
        if not read1:
            break
        
        if ".fastq" in read1:
            sampleName = ((read1.split("/"))[-1].split("_1.fastq"))[0]
        if ".fq" in read1:
            sampleName = ((read1.split("/"))[-1].split("_1.fq"))[0]
        
        if not ".fastq" in read1 and not ".fq" in read1:
            print "Reads file must have the suffix .fastq or .fq. Now exit...."
            exit(0)

        print "Analising sample",sampleName
        print "Aigning reads to reference"
        os.system("bowtie2 --end-to-end -1 "+read1+" -2 "+read2+" -x reference -S alignment.sam -p "+numThreads)
        print "Converting sam to bam"
        os.system("samtools view -bS -h alignment.sam >alignment.bam")
        print "Sorting bam"
        os.system("samtools sort -o "+sampleName+"_alignment_sorted.bam alignment.bam")
        print "Calling snps with lofreq"
        #Analyze snps with lowfreq
        os.system("lofreq call -f "+referenceFile+" -o "+sampleName+".vcf "+sampleName+"_alignment_sorted.bam")
        
        #Analyze indes using the GATK pipeline
        os.system("samtools view -bF 4 "+sampleName+"_alignment_sorted.bam >mapped.bam")
        os.system("~/Software/jre1.8.0_191/bin/java -jar -XX:ParallelGCThreads=16  ~/Software/picard.jar AddOrReplaceReadGroups I=mapped.bam O=rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=Ilumina RGPU=machine RGSM=Consensus")
        os.system("~/Software/jre1.8.0_191/bin/java -jar -XX:ParallelGCThreads=16  ~/Software/picard.jar MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics")
        os.system("~/Software/jre1.8.0_191/bin/java -jar ~/Software/picard.jar CreateSequenceDictionary R="+referenceFile)
        os.system("samtools faidx "+referenceFile)
        os.system("~/Software/jre1.8.0_191/bin/java -jar  ~/Software/gatk3.8/GenomeAnalysisTK.jar -T  HaplotypeCaller -R "+referenceFile+" -I dedupped.bam  -o output.vcf -A StrandAlleleCountsBySample")
        os.system("mv output.vcf "+sampleName+"_indels.vcf")
               
        
        os.system("mv "+sampleName+"_alignment_sorted.bam "+outputFolder+"/")
        os.system("mv *.vcf "+outputFolder+"/")
        os.system("rm alignment* -f")
        


