from Bio import SeqIO
import biomodule as bm
import sys
import time
import os
from os import listdir
from os.path import isfile, join

genomeName = sys.argv[1]
outputFolder = sys.argv[2]
prot2map = []

suffixName = (genomeName.split("/"))[-1]



gffFile = open(suffixName+"_annotation.gff","w") #To change from the command line
warnFile = open(suffixName+"_annotationWarnings.txt","w") #To change from the command line
cdsFile = open(suffixName+"_cds.fasta","w") #To change from the command line
protFile = open(suffixName+"_proteins.fasta","w") #To change from the command line


os.system("rm -f ./proteinDB/._*")


def checkCDS(cds):
    for a in range(0,len(cds)-3,+3):
        if cds[a:a+3] == "TAA" or cds[a:a+3] == "TGA" or cds[a:a+3] == "TAG":
            return "Stop codon in CDS"
    if not float(len(cds))%3.0==0:
        return "Not multiple of 3 coding sequence length"



for seq_record in SeqIO.parse(genomeName,"fasta"):
    genomeSeq = str(seq_record.seq)
    assemblyName = str(seq_record.id)


onlyfiles = [f for f in listdir("./proteinDB/") if isfile(join("./proteinDB/", f))]
for f in onlyfiles:
    if f[-13:] == "_models.fasta":
        prot2map.append(f)



for f in prot2map:
    numCodonRefines = 0
    #Find the best model for the gene *********************************************************
    locus = f.replace("_models.fasta","")
    protSeqs = SeqIO.to_dict(SeqIO.parse("./proteinDB/"+f,"fasta"))
    print "Choosing best match for protein",f

    os.system("blastx -query "+genomeName+" -db ./proteinDB/"+f+" -outfmt 6 | sort   -k 12rn,12rn > outputBlast.txt")
    blastFile = open("outputBlast.txt")
    bestProt = ((blastFile.readline()).split("\t"))[1]
    print "Best match in the following protein:",bestProt

    sequence = str(protSeqs[bestProt].seq)
    tempFasta = open("tempFasta.fasta","w")
    tempFasta.write(">"+locus+"\n"+str(protSeqs[bestProt].seq)+"\n")
    tempFasta.close()

    #Run Exonerate on the best model *********************************************************
    os.system("exonerate --model protein2genome tempFasta.fasta "+genomeName+" --showtargetgff -s 0 -n 1 --forcegtag  >outputExonerate")
   
   
   
    #Check Exonerate output ***************************************************************** 
    #Check the proteins gave a match in the target genome
    exResult = open("outputExonerate")
    line = exResult.readline().rstrip()
    while not "Query range:" in line:
        line = exResult.readline().rstrip()
        if line is None:
            print "WARNING no protein found for ", locus
            warnings.append("Missing gene: locus "+locus+" did not provide any alignment")
    exResult.close()

    
    #Reconstruct Exons  **********************************************************************
    exResult = open("outputExonerate")
    while not line == "# --- START OF GFF DUMP ---":
        line = exResult.readline().rstrip()
        if line is None:
            #print "WARNING no protein found for ", locus
            warnings.append("Missing gene: locus "+locus+" did not provide any alignment")
    for a in range(10):
        line = exResult.readline().rstrip()

    gene = {}
    exon = {}
    warnings = []
    #Collect the exonerate output
    while not line == "# --- END OF GFF DUMP ---":
        line = exResult.readline().rstrip()
        if line is None:
            #print "WARNING no protein found for ", locus
            warnings.append("Missing gene: locus "+locus+" did not provide any alignment")
        fields = line.split("\t")
        if not line == "# --- END OF GFF DUMP ---":
            if fields[2]=="gene":
                if not locus in gene:
                    gene[locus] = (fields[3],fields[4],fields[6])
                else:
                    if (int(fields[4]) - int(fields[3])) > (int(int(gene[locus][1])) - int(int(gene[locus][0]))):
                        gene[locus] = (fields[3],fields[4],fields[6])
            
            if fields[2] == "exon":
                if not locus in exon:
                    exon[locus] = []
                exon[locus].append((fields[3],fields[4],fields[6]))
                if "frameshifts" in fields[8]:
                    warnings.append("Frameshifts in exon "+str(fields[3])+" "+str(fields[4])+" "+fields[8])
   
    newList = sorted(exon[locus], key=lambda x: x[1])
    exon[locus] = newList


    #Reconstruct CDS  ***********************************************************
    cdsSeq = "" 
    if exon[locus][0][2]=="+":   #************* Positive strand
        for item in exon[locus]:
            cdsSeq+=genomeSeq[int(item[0])-1:int(item[1])]
        cdsSeq += genomeSeq[int(item[1]):int(item[1])+3]

    else: # *********************** Negative Strand
        for a in range(len(exon[locus])-1,-1,-1):
            cdsSeq+=bm.reverseComplement(genomeSeq[ int(exon[locus][a][0])-1: int(exon[locus][a][1]) ])
        cdsSeq += bm.reverseComplement(genomeSeq[int(exon[locus][a][0])-4 :int(exon[locus][a][0])-1])


    notes = ""

    # Check CDS integrity *********************************************************
    foundStartCodon = True
    foundStopCodon = True
    if not  cdsSeq[:3]=="ATG" or not (cdsSeq[-3:]=="TGA" or cdsSeq[-3:]=="TAA" or cdsSeq[-3:]=="TAG" ):

        notes = "\n"+locus+"\n"
        



        # Look for ATG at the beginning of the sequence or closely ********************
        if exon[locus][0][2]=="+":   #************* Positive strand
            if not cdsSeq[:3]=="ATG":
                foundStartCodon = False
                for a in range(len(sequence)-len(cdsSeq)/3+30):
                    newStart = genomeSeq[int(exon[locus][0][0])-1-a*3-3:int(exon[locus][0][0])-1-a*3]
                    if newStart == "ATG":
                        exon[locus][0]=(int(exon[locus][0][0])-a*3-3,exon[locus][0][1],exon[locus][0][2])
                        gene[locus] = (int(exon[locus][0][0])-a*3-3,int(gene[locus][1]),gene[locus][2])
                        cdsSeq = ""
                        for item in exon[locus]:
                            cdsSeq+=genomeSeq[int(item[0])-1:int(item[1])]     
                        foundStartCodon = True
                        notes += "- Start codon refined  "+str(a)+" codons upstream\n"
                        numCodonRefines = a
                        break
                    if newStart == "TGA" or newStart == "TAA" or newStart=="TAG":
                        foundStartCodon = False
                        notes += "- No Start codon found "
                        break
                #If the new start codon was not found in the region upstream then the downstream region is searched
                for a in range(len(sequence)-len(cdsSeq)/3+30):
                    newStart = genomeSeq[int(exon[locus][0][0])-1+a*3+3:int(exon[locus][0][0])-1+a*3]
                    if newStart == "ATG":
                        exon[locus][0]=(int(exon[locus][0][0])+a*3+3,exon[locus][0][1],exon[locus][0][2])
                        gene[locus] = (int(exon[locus][0][0])+a*3+3,int(gene[locus][1]),gene[locus][2])
                        cdsSeq = ""
                        for item in exon[locus]:
                            cdsSeq+=genomeSeq[int(item[0])-1:int(item[1])]     
                        foundStartCodon = True
                        notes += "- Start codon refined  "+str(a)+" codons upstream\n"
                        numCodonRefines = a
                        break
                    if newStart == "TGA" or newStart == "TAA" or newStart=="TAG":
                        foundStartCodon = False
                        notes += "- No Start codon found "
                        break

        else: # *********************** Negative Strand
            if not cdsSeq[:3]=="ATG":
                foundStartCodon = False
                for a in range(len(sequence)-len(cdsSeq)/3+30):
                    #print "New start codons"
                    newStart = bm.reverseComplement(genomeSeq[int(exon[locus][-1][1])+a*3:int(exon[locus][-1][1])+a*3+3])
                    #print newStart
                    if newStart == "ATG":
                        exon[locus][-1]=(int(exon[locus][-1][0]),int(exon[locus][-1][1])+a*3+3,exon[locus][-1][2])
                        gene[locus] = (int(gene[locus][0]),int(exon[locus][-1][1])+a*3+3,gene[locus][2])
                        cdsSeq = ""
                        for a in range(len(exon[locus])-1,-1,-1):
                            cdsSeq+=bm.reverseComplement(genomeSeq[ int(exon[locus][a][0])-1: int(exon[locus][a][1]) ])   
                        foundStartCodon = 1
                        notes += "- Start codon refined  "+str(a)+" codons upstream\n"
                        numCodonRefines = a
                        break
                    if newStart == "TGA" or newStart == "TAA" or newStart=="TAG":
                        foundStartCodon = False
                        notes += "- No Start codon found "
                        break
            #If the new start codon was not found in the region upstream then the downstream region is searched
                for a in range(len(sequence)-len(cdsSeq)/3+30):
                    #print "New start codons"
                    newStart = bm.reverseComplement(genomeSeq[int(exon[locus][-1][1])-a*3:int(exon[locus][-1][1])-a*3-3])
                    #print newStart
                    if newStart == "ATG":
                        exon[locus][-1]=(int(exon[locus][-1][0]),int(exon[locus][-1][1])-a*3-3,exon[locus][-1][2])
                        gene[locus] = (int(gene[locus][0]),int(exon[locus][-1][1])-a*3-3,gene[locus][2])
                        cdsSeq = ""
                        for a in range(len(exon[locus])-1,-1,-1):
                            cdsSeq+=bm.reverseComplement(genomeSeq[ int(exon[locus][a][0])-1: int(exon[locus][a][1]) ])   
                        foundStartCodon = 1
                        notes += "- Start codon refined  "+str(a)+" codons upstream\n"
                        numCodonRefines = a
                        break
                    if newStart == "TGA" or newStart == "TAA" or newStart=="TAG":
                        foundStartCodon = False
                        notes += "- No Start codon found "
                        break

            


        # Look for Stop codon at the end of the sequence or closely ********************
        if exon[locus][0][2]=="+":   #************* Positive strand
            if not cdsSeq[-3:]=="TGA" and not cdsSeq[-3:]=="TAA" and not cdsSeq[-3:]=="TAG":
                foundStopCodon = False
                for a in range(len(sequence)-len(cdsSeq)/3+30):
                    newStop = genomeSeq[int(exon[locus][-1][1])+a*3:int(exon[locus][-1][1])+a*3+3]
                    if newStop == "TAA" or newStop=="TGA" or newStop=="TAG":
                        exon[locus][0]=(int(exon[locus][0][0]),int(exon[locus][0][1])+a*3+3,exon[locus][0][2])
                        gene[locus] = (int(gene[locus][0]), int(exon[locus][0][1])+a*3+3, gene[locus][2])
                        cdsSeq = ""
                        for item in exon[locus]:
                            cdsSeq+=genomeSeq[int(item[0])-1:int(item[1])]     
                        foundStopCodon = True
                        notes += "- Stop codon refined "+str(a)+" codon upstream\n"
                        break


        else: # *********************** Negative Strand
            if not cdsSeq[-3:]=="TGA" and not cdsSeq[-3:]=="TAA" and not cdsSeq[-3:]=="TAG":
                foundStopCodon = False
                for a in range(len(sequence)-len(cdsSeq)/3+30):
                    #print "New Stop codons"
                    newStop = bm.reverseComplement(genomeSeq[int(exon[locus][0][0])-1-a*3-3:int(exon[locus][0][0])-1-a*3])
                    #print newStop
                    if newStop == "TAA" or newStop=="TGA" or newStop=="TAG":
                        exon[locus][0]=(int(exon[locus][0][0])-a*3-3,exon[locus][0][1],exon[locus][0][2])
                        gene[locus] =(int(exon[locus][0][0])-a*3-3, int(gene[locus][1]),gene[locus][2])
                        cdsSeq = ""
                        for b in range(len(exon[locus])-1,-1,-1):
                            cdsSeq+=bm.reverseComplement(genomeSeq[ int(exon[locus][b][0])-1: int(exon[locus][b][1]) ])    
                        foundStopCodon = True
                        notes += "- Stop codon refined "+str(a)+" codon upstream\n"
                        break

        warnFile.write(notes+"\n")
    

    if foundStartCodon == True and foundStopCodon == True and numCodonRefines <5: # Write gff and cds file ********************
        if  notes == "":
            notes = "\n"+locus+"\n"

        #  ******************* Check CDS integrity
        cdsGood = True
        for a in range(0,len(cdsSeq)-3,+3):
            if cdsSeq[a:a+3] == "TAA" or cdsSeq[a:a+3] == "TGA" or cdsSeq[a:a+3] == "TAG":
                # Check if the shorter sequence is compatible with one of the models
                newProtLen = float(a)/3.0
                for protein in protSeqs:
                    if newProtLen / float(len(protSeqs[protein].seq)) >= 0.9:
                        if exon[locus][0][2]=="+":  #Check it if the strand is positive
                            newmRNALength = 0
                            newExonSet = {}
                            if not locus in newExonSet:
                                newExonSet[locus] = []
                            for item in exon[locus]:
                                if int(item[1])-int(item[0]) + newmRNALength > a+3:
                                    newExonSet[locus].append((int(item[0]),int(item[0]) + a+3 - newmRNALength -1 ,item[2]))
                                    exon[locus] = newExonSet[locus]
                                    gene[locus] = (int(gene[locus][0]),int(item[0]) + a+3 - newmRNALength -1 ,gene[locus][2])
                                    break
                                else:
                                    newmRNALength += int(item[1]) - int(item[0])
                                    newExonSet[locus].append((int(item[0]),int(item[1]),item[2]))
                            cdsSeq = cdsSeq[:a+3]
                        else:  #Check it if the strand is negative
                            newmRNALength = 0
                            newExonSet = {}
                            if not locus in newExonSet:
                                newExonSet[locus] = []
                            for a in range(len(exon[locus])-1,-1,-1):
                                if int(item[1])-int(item[0]) + newmRNALength > a+3:
                                    newExonSet[locus].append((int(item[1]) - a -3 + newmRNALength, int(item[1]),item[2]))
                                    exon[locus] = newExonSet[locus]
                                    gene[locus] = (int(item[1]) - a -3 + newmRNALength, int(gene[locus][1]), gene[locus][2])
                                    break

                                else:
                                    newmRNALength = int(item[1]) - int(item[0])
                                    newExonSet[locus].append((int(item[0]),int(item[1]),item[2]))
                            cdsSeq = cdsSeq[:a+3]

                        break
                    else:
                        cdsGood = False
                        gffNote = "note=Stop codon interrupts coding sequence. "
                        notes += "The coding sequence is interrupted by a stop codon\n"
                        break
            
        if not len(cdsSeq)%3 == 0:
            cdsGood = False
            notes += "The coding sequence is interrupted by a stop codon\n"
            if not gffNote == "" :
                gffNote = "note= insertions/deletions lead to broken codons."
            else:
                gffNote +="Insertions/deletions lead to broken codons."
            

            

        if cdsGood == True:  # CDS passed quality check
            if exon[locus][0][2]=="+":   #************* Positive strand
                gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(int(gene[locus][0]))+"\t"+str(int(gene[locus][1])+3)+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+"\n")
                gffFile.write(assemblyName+"\texonerate\t"+"mRNA"+"\t"+str(int(gene[locus][0]))+"\t"+str(int(gene[locus][1])+3)+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+"\n")
                numExon = 1
                for item in exon[locus]:
                    gffFile.write(assemblyName+"\texonerate\t"+"CDS"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+"\n")
                    numExon += 1
                

            else: # *********************** Negative Strand
                gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(int(gene[locus][0])-3)+"\t"+str(int(gene[locus][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+"\n")
                gffFile.write(assemblyName+"\texonerate\t"+"mRNA"+"\t"+str(int(gene[locus][0])-3)+"\t"+str(int(gene[locus][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+"\n")
                numExon = 1
                for item in exon[locus]:
                    gffFile.write(assemblyName+"\texonerate\t"+"CDS"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+"\n")
                    numExon += 1

            cdsFile.write(">"+locus+" +\n"+cdsSeq+"\n")

        else: # CDS DID NOT passed quality check
            warnFile.write(notes+"\n")
            gffNote += "Pseudo "
            if exon[locus][0][2]=="+":   #************* Positive strand
                gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(int(gene[locus][0]))+"\t"+str(int(gene[locus][1])+3)+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+"\t"+gffNote+"-gene\n")
                gffFile.write(assemblyName+"\texonerate\t"+"misc_feature"+"\t"+str(int(gene[locus][0]))+"\t"+str(int(gene[locus][1])+3)+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+";"+gffNote+"-mRNA\n")
                numExon = 1
                for item in exon[locus]:
                    gffFile.write(assemblyName+"\texonerate\t"+"misc_feature"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+";"+gffNote+"-CDS\n")
                    numExon += 1
                

            else: # *********************** Negative Strand
                gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(int(gene[locus][0])-3)+"\t"+str(int(gene[locus][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+";"+gffNote+"-gene\n")
                gffFile.write(assemblyName+"\texonerate\t"+"misc_feature"+"\t"+str(int(gene[locus][0])-3)+"\t"+str(int(gene[locus][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+";"+gffNote+"-mRNA\n")
                numExon = 1
                for item in exon[locus]:
                    gffFile.write(assemblyName+"\texonerate\t"+"misc_feature"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+";"+gffNote+"-CDS\n")
                    numExon += 1

            cdsFile.write(">"+locus+" +\n"+cdsSeq+"\n")




    

    if foundStartCodon == False or foundStopCodon == False or numCodonRefines >=5:
        print "Annotation needs refinement"
        exResult.close()
        #  ***************************************************************************
        #  ************************* Annotation refinement ***************************
        #  ***************************************************************************
        

        os.system("exonerate --model protein2genome tempFasta.fasta "+genomeName+" --showtargetgff -s 0 -n 1 --refine full --forcegtag >outputExonerate")
        
        #Check Exonerate output ***************************************************************** 
        #Check the proteins gave a match in the target genome
        exResult = open("outputExonerate")
        line = exResult.readline().rstrip()
        while not "Query range:" in line:
            line = exResult.readline().rstrip()
            if line is None:
                print "WARNING no protein found for ", locus
                warnings.append("Missing gene: locus "+locus+" did not provide any alignment")
        exResult.close()

        
        #Reconstruct Exons  **********************************************************************
        exResult = open("outputExonerate")
        while not line == "# --- START OF GFF DUMP ---":
            line = exResult.readline().rstrip()
            if line is None:
                #print "WARNING no protein found for ", locus
                warnings.append("Missing gene: locus "+locus+" did not provide any alignment")
        for a in range(10):
            line = exResult.readline().rstrip()

        gene = {}
        exon = {}
        warnings = []
        #Collect the exonerate output
        while not line == "# --- END OF GFF DUMP ---":
            line = exResult.readline().rstrip()
            if line is None:
                #print "WARNING no protein found for ", locus
                warnings.append("Missing gene: locus "+locus+" did not provide any alignment")
            fields = line.split("\t")
            if not line == "# --- END OF GFF DUMP ---":
                if fields[2]=="gene":
                    if not locus in gene:
                        gene[locus] = (fields[3],fields[4],fields[6])
                    else:
                        if (int(fields[4]) - int(fields[3])) > (int(int(gene[locus][1])) - int(int(gene[locus][0]))):
                            gene[locus] = (fields[3],fields[4],fields[6])
                
                if fields[2] == "exon":
                    if not locus in exon:
                        exon[locus] = []
                    exon[locus].append((fields[3],fields[4],fields[6]))
                    if "frameshifts" in fields[8]:
                        warnings.append("Frameshifts in exon "+str(fields[3])+" "+str(fields[4])+" "+fields[8])
    
        newList = sorted(exon[locus], key=lambda x: x[1])
        exon[locus] = newList


        #Reconstruct CDS  ***********************************************************
        cdsSeq = "" 
        if exon[locus][0][2]=="+":   #************* Positive strand
            for item in exon[locus]:
                cdsSeq+=genomeSeq[int(item[0])-1:int(item[1])]
            cdsSeq += genomeSeq[int(item[1]):int(item[1])+3]

        else: # *********************** Negative Strand
            for a in range(len(exon[locus])-1,-1,-1):
                cdsSeq+=bm.reverseComplement(genomeSeq[ int(exon[locus][a][0])-1: int(exon[locus][a][1]) ])
            cdsSeq += bm.reverseComplement(genomeSeq[int(exon[locus][a][0])-4 :int(exon[locus][a][0])-1])




        # Check CDS completness *********************************************************
        foundStartCodon = True
        foundStopCodon = True
        if not  cdsSeq[:3]=="ATG" or not (cdsSeq[-3:]=="TGA" or cdsSeq[-3:]=="TAA" or cdsSeq[-3:]=="TAG" ):

            notes = locus+" - after refinement \n"
            



            # Look for ATG at the beginning of the sequence or closely ********************
            if exon[locus][0][2]=="+":   #************* Positive strand
                if not cdsSeq[:3]=="ATG":
                    foundStartCodon = False
                    for a in range(len(sequence)-len(cdsSeq)/3+30):
                        newStart = genomeSeq[int(exon[locus][0][0])-1-a*3-3:int(exon[locus][0][0])-1-a*3]
                        if newStart == "ATG":
                            exon[locus][0]=(int(exon[locus][0][0])-a*3-3,exon[locus][0][1],exon[locus][0][2])
                            gene[locus] = (int(exon[locus][0][0])-a*3-3,int(gene[locus][1]),gene[locus][2])
                            cdsSeq = ""
                            for item in exon[locus]:
                                cdsSeq+=genomeSeq[int(item[0])-1:int(item[1])]     
                            foundStartCodon = True
                            notes += "- Start codon refined  "+str(a)+" codons upstream\n"
                            numCodonRefines = a
                            break
                        if newStart == "TGA" or newStart == "TAA" or newStart=="TAG":
                            foundStartCodon = False
                            notes += "- No Start codon found "
                            break
            #If the new start codon was not found in the region upstream then the downstream region is searched
                    for a in range(len(sequence)-len(cdsSeq)/3+30):
                        newStart = genomeSeq[int(exon[locus][0][0])-1+a*3+3:int(exon[locus][0][0])-1+a*3]
                        if newStart == "ATG":
                            exon[locus][0]=(int(exon[locus][0][0])+a*3+3,exon[locus][0][1],exon[locus][0][2])
                            gene[locus] = (int(exon[locus][0][0])+a*3+3,int(gene[locus][1]),gene[locus][2])
                            cdsSeq = ""
                            for item in exon[locus]:
                                cdsSeq+=genomeSeq[int(item[0])-1:int(item[1])]     
                            foundStartCodon = True
                            notes += "- Start codon refined  "+str(a)+" codons upstream\n"
                            numCodonRefines = a
                            break
                        if newStart == "TGA" or newStart == "TAA" or newStart=="TAG":
                            foundStartCodon = False
                            notes += "- No Start codon found "
                            break

            else: # *********************** Negative Strand
                if not cdsSeq[:3]=="ATG":
                    foundStartCodon = False
                    for a in range(len(sequence)-len(cdsSeq)/3+30):
                        #print "New start codons"
                        newStart = bm.reverseComplement(genomeSeq[int(exon[locus][-1][1])+a*3:int(exon[locus][-1][1])+a*3+3])
                        #print newStart
                        if newStart == "ATG":
                            exon[locus][-1]=(int(exon[locus][-1][0]),int(exon[locus][-1][1])+a*3+3,exon[locus][-1][2])
                            gene[locus] = (int(gene[locus][0]),int(exon[locus][-1][1])+a*3+3,gene[locus][2])
                            cdsSeq = ""
                            for a in range(len(exon[locus])-1,-1,-1):
                                cdsSeq+=bm.reverseComplement(genomeSeq[ int(exon[locus][a][0])-1: int(exon[locus][a][1]) ])   
                            foundStartCodon = 1
                            notes += "- Start codon refined  "+str(a)+" codons upstream\n"
                            numCodonRefines = a
                            break
                        if newStart == "TGA" or newStart == "TAA" or newStart=="TAG":
                            foundStartCodon = False
                            notes += "- No Start codon found "
                            break
                #If the new start codon was not found in the region upstream then the downstream region is searched
                    for a in range(len(sequence)-len(cdsSeq)/3+30):
                        #print "New start codons"
                        newStart = bm.reverseComplement(genomeSeq[int(exon[locus][-1][1])-a*3-3:int(exon[locus][-1][1])-a*3])
                        #print newStart



                        if newStart == "ATG":
                            exon[locus][-1]=(int(exon[locus][-1][0]),int(exon[locus][-1][1])-a*3-3,exon[locus][-1][2])
                            gene[locus] = (int(gene[locus][0]),int(exon[locus][-1][1])-a*3-3,gene[locus][2])
                            cdsSeq = ""
                            for a in range(len(exon[locus])-1,-1,-1):
                                cdsSeq+=bm.reverseComplement(genomeSeq[ int(exon[locus][a][0])-1: int(exon[locus][a][1]) ])   
                            foundStartCodon = 1
                            notes += "- Start codon refined  "+str(a)+" codons upstream\n"
                            numCodonRefines = a
                            break
                        if newStart == "TGA" or newStart == "TAA" or newStart=="TAG":
                            foundStartCodon = False
                            notes += "- No Start codon found "
                            break

            # Look for Stop codon at the end of the sequence or closely ********************
            if exon[locus][0][2]=="+":   #************* Positive strand
                if not cdsSeq[-3:]=="TGA" and not cdsSeq[-3:]=="TAA" and not cdsSeq[-3:]=="TAG":
                    foundStopCodon = False
                    for a in range(len(sequence)-len(cdsSeq)/3+30):
                        newStop = genomeSeq[int(exon[locus][-1][1])+a*3:int(exon[locus][-1][1])+a*3+3]
                        if newStop == "TAA" or newStop=="TGA" or newStop=="TAG":
                            exon[locus][0]=(int(exon[locus][0][0]),int(exon[locus][0][1])+a*3+3,exon[locus][0][2])
                            gene[locus] = (int(gene[locus][0]), int(exon[locus][0][1])+a*3+3, gene[locus][2])
                            cdsSeq = ""
                            for item in exon[locus]:
                                cdsSeq+=genomeSeq[int(item[0])-1:int(item[1])]     
                            foundStopCodon = True
                            notes += "- Stop codon refined "+str(a)+" codon upstream\n"
                            break


            else: # *********************** Negative Strand
                if not cdsSeq[-3:]=="TGA" and not cdsSeq[-3:]=="TAA" and not cdsSeq[-3:]=="TAG":
                    foundStopCodon = 0
                    for a in range(len(sequence)-len(cdsSeq)/3+30):
                        #print "New Stop codons"
                        newStop = bm.reverseComplement(genomeSeq[int(exon[locus][0][0])-1-a*3-3:int(exon[locus][0][0])-1-a*3])
                        #print newStop
                        if newStop == "TAA" or newStop=="TGA" or newStop=="TAG":
                            exon[locus][0]=(int(exon[locus][0][0])-a*3-3,exon[locus][0][1],exon[locus][0][2])
                            gene[locus] =(int(exon[locus][0][0])-a*3-3, int(gene[locus][1]),gene[locus][2])
                            cdsSeq = ""
                            for b in range(len(exon[locus])-1,-1,-1):
                                cdsSeq+=bm.reverseComplement(genomeSeq[ int(exon[locus][b][0])-1: int(exon[locus][b][1]) ])    
                            foundStopCodon = 1
                            notes += "- Stop codon refined "+str(a)+" codon upstream\n"
                            break

            warnFile.write(notes+"\n")



        if foundStartCodon == True and foundStopCodon == True: # Write gff and cds file ********************
            cdsGood = True
           #  ******************* Check CDS integrity
            for a in range(0,len(cdsSeq)-3,+3):
                if cdsSeq[a:a+3] == "TAA" or cdsSeq[a:a+3] == "TGA" or cdsSeq[a:a+3] == "TAG":
                    # Check if the shorter sequence is compatible with one of the models
                    newProtLen = float(a)/3.0
                    for protein in protSeqs:
                        if newProtLen / float(len(protSeqs[protein].seq)) >= 0.9:
                            if exon[locus][0][2]=="+":  #Check it if the strand is positive
                                newmRNALength = 0
                                newExonSet = {}
                                if not locus in newExonSet:
                                    newExonSet[locus] = []
                                for item in exon[locus]:
                                    if int(item[1])-int(item[0]) + newmRNALength > a+3:
                                        newExonSet[locus].append((int(item[0]),int(item[0]) + a+3 - newmRNALength -1 ,item[2]))
                                        exon[locus] = newExonSet[locus]
                                        gene[locus] = (int(gene[locus][0]),int(item[0]) + a+3 - newmRNALength -1 ,gene[locus][2])
                                        break
                                    else:
                                        newmRNALength += int(item[1]) - int(item[0])
                                        newExonSet[locus].append((int(item[0]),int(item[1]),item[2]))
                                cdsSeq = cdsSeq[:a+3]
                            else:  #Check it if the strand is negative
                                newmRNALength = 0
                                newExonSet = {}
                                if not locus in newExonSet:
                                    newExonSet[locus] = []
                                for a in range(len(exon[locus])-1,-1,-1):
                                    if int(item[1])-int(item[0]) + newmRNALength > a+3:
                                        newExonSet[locus].append((int(item[1]) - a -3 + newmRNALength, int(item[1]),item[2]))
                                        exon[locus] = newExonSet[locus]
                                        gene[locus] = (int(item[1]) - a -3 + newmRNALength, int(gene[locus][1]), gene[locus][2])
                                        break

                                    else:
                                        newmRNALength = int(item[1]) - int(item[0])
                                        newExonSet[locus].append((int(item[0]),int(item[1]),item[2]))
                                cdsSeq = cdsSeq[:a+3]

                            break
                        else:
                            cdsGood = False
                            gffNote = "note=Stop codon interrupts coding sequence. "
                            notes += "The coding sequence is interrupted by a stop codon\n"
                            break

            if not len(cdsSeq)%3 == 0:
                cdsGood = False
                notes += "Insertions/deletions lead to broken codons\n"
                if not gffNote == "" :
                    gffNote = "note= insertions/deletions lead to broken codons."
                else:
                    gffNote +="Insertions/deletions lead to broken codons."
                

                

            if cdsGood == True:  # CDS passed quality check
                notes += "\n"+ locus+" - after refinement \n"+"Good CDS!"
                warnFile.write(notes+"\n")
                if exon[locus][0][2]=="+":   #************* Positive strand
                    gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(int(gene[locus][0]))+"\t"+str(int(gene[locus][1])+3)+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+"\n")
                    gffFile.write(assemblyName+"\texonerate\t"+"mRNA"+"\t"+str(int(gene[locus][0]))+"\t"+str(int(gene[locus][1])+3)+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+"\n")
                    numExon = 1
                    for item in exon[locus]:
                        gffFile.write(assemblyName+"\texonerate\t"+"CDS"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+"\n")
                        numExon += 1
                    

                else: # *********************** Negative Strand
                    gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(int(gene[locus][0])-3)+"\t"+str(int(gene[locus][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+"\n")
                    gffFile.write(assemblyName+"\texonerate\t"+"mRNA"+"\t"+str(int(gene[locus][0])-3)+"\t"+str(int(gene[locus][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+"\n")
                    numExon = 1
                    for item in exon[locus]:
                        gffFile.write(assemblyName+"\texonerate\t"+"CDS"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+"\n")
                        numExon += 1

                cdsFile.write(">"+locus+" +\n"+cdsSeq+"\n")

            else: # CDS DID NOT passed quality check
                warnFile.write(notes+"\n")
                gffNote += "Pseudo "
                if exon[locus][0][2]=="+":   #************* Positive strand
                    gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(int(gene[locus][0]))+"\t"+str(int(gene[locus][1])+3)+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+";"+gffNote+"-gene\n")
                    gffFile.write(assemblyName+"\texonerate\t"+"misc_feature"+"\t"+str(int(gene[locus][0]))+"\t"+str(int(gene[locus][1])+3)+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+";"+gffNote+"-mRNA\n")
                    numExon = 1
                    for item in exon[locus]:
                        gffFile.write(assemblyName+"\texonerate\t"+"misc_feature"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+";"+gffNote+"-CDS\n")
                        numExon += 1
                    

                else: # *********************** Negative Strand
                    gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(int(gene[locus][0])-3)+"\t"+str(int(gene[locus][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+";"+gffNote+"-gene\n")
                    gffFile.write(assemblyName+"\texonerate\t"+"misc_feature"+"\t"+str(int(gene[locus][0])-3)+"\t"+str(int(gene[locus][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+";"+gffNote+"-mRNA\n")
                    numExon = 1
                    for item in exon[locus]:
                        gffFile.write(assemblyName+"\texonerate\t"+"misc_feature"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+";"+gffNote+"-CDS\n")
                        numExon += 1

                cdsFile.write(">"+locus+" " +gffNote+"-CDS\n"+cdsSeq+"\n")

        else:
            gffNote = "Possible pseudogene"
            if exon[locus][0][2]=="+":   #************* Positive strand
                gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(int(gene[locus][0]))+"\t"+str(int(gene[locus][1])+3)+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+";"+gffNote+"-gene\n")
                gffFile.write(assemblyName+"\texonerate\t"+"misc_feature"+"\t"+str(int(gene[locus][0]))+"\t"+str(int(gene[locus][1])+3)+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+";"+gffNote+"-mRNA\n")
                numExon = 1
                for item in exon[locus]:
                    gffFile.write(assemblyName+"\texonerate\t"+"misc_feature"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+";"+gffNote+"-CDS\n")
                    numExon += 1
                    

            else: # *********************** Negative Strand
                gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(int(gene[locus][0])-3)+"\t"+str(int(gene[locus][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+";"+gffNote+"-gene\n")
                gffFile.write(assemblyName+"\texonerate\t"+"misc_feature"+"\t"+str(int(gene[locus][0])-3)+"\t"+str(int(gene[locus][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+";"+gffNote+"-mRNA\n")
                numExon = 1
                for item in exon[locus]:
                    gffFile.write(assemblyName+"\texonerate\t"+"misc_feature"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+";"+gffNote+"-CDS\n")
                    numExon += 1

            cdsFile.write(">"+locus+" " +gffNote+"-CDS\n"+cdsSeq+"\n")
           



os.system("rm tempFasta.fasta outputExonerate outputBlast.txt -f")   
cdsFile.close()
gffFile.close()
warnFile.close()

#Translate valid cds in proteins 

sequences = SeqIO.to_dict(SeqIO.parse(suffixName+"_cds.fasta","fasta"))

for seq in sequences:
    if not "pseudo" in str(sequences[seq].description):
        protSeq = (sequences[seq].seq).translate()
        protFile.write(">"+str(sequences[seq].description)+"\n"+str(protSeq)+"\n")

protFile.close()


os.system("mv "+suffixName+"* "+outputFolder)


        