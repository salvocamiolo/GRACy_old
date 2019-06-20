try:
    import Tkinter as tk
except ImportError:
    import tkinter as tk

try:
    import ttk
    py3 = False
except ImportError:
    import tkinter.ttk as ttk
    py3 = True

import Tkinter, Tkconstants, tkFileDialog
import sys
from Bio import SeqIO
import biomodule as bm

import time
import os
from os import listdir
from os.path import isfile, join

installationDirectory = sys.argv[1]

def vp_start_gui():
    '''Starting point when module is the main routine.'''
    global val, w, root

    root = tk.Tk()

   

    top = Toplevel1 (root)
    root.mainloop()



def create_Toplevel1(root, *args, **kwargs):
    
    '''Starting point when module is imported by another program.'''
    global w, w_win, rt
    rt = root
    w = tk.Toplevel (root)
    VirHosFilt_support.set_Tk_var()
    top = Toplevel1 (w)
    VirHosFilt_support.init(w, top, *args, **kwargs)
    return (w, top)

def destroy_Toplevel1():
    global w
    w.destroy()
    w = None

class Toplevel1:
    def __init__(self, top=None):
        
        def openInputFile():
            inputFile = tkFileDialog.askopenfilename(initialdir = "./",title = "Select an input file")
            self.inputFileEntry.delete(0,tk.END)
            self.inputFileEntry.insert(0,inputFile)

        def openOutputFolder():
            outputFolder = tkFileDialog.askdirectory(initialdir = "./",title = "Select folder")
            self.outputFolderEntry.delete(0,tk.END)
            self.outputFolderEntry.insert(0,outputFolder)

        def exitProgram():
            exit()

        def checkCDS(cds):
            for a in range(0,len(cds)-3,+3):
                if cds[a:a+3] == "TAA" or cds[a:a+3] == "TGA" or cds[a:a+3] == "TAG":
                    return "Stop codon in CDS"
            if not float(len(cds))%3.0==0:
                return "Not multiple of 3 coding sequence length"

        def mainAnnotationAlgorithm():
            if self.inputFileEntry.get() == "Please select a genomes list...." or self.outputFolderEntry.get() == "Please select an output folder....":
                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "Please select a valid input file and output folder....\n")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()
            else:
                #**************************************************************************
                #***************** Main Annotation Algorithm Start ************************
                #**************************************************************************
                infile = open(self.inputFileEntry.get())
                outputFolder = self.outputFolderEntry.get()
                while True:
                    file2Annotate = infile.readline().rstrip()
                    if not file2Annotate:
                        break
                    file2AnnotateName = (file2Annotate.split("/"))[-1]
                    self.logArea.configure(state='normal')
                    self.logArea.insert(tk.END, "Starting annotation on sample "+file2AnnotateName+"....\n")
                    self.logArea.see(tk.END)
                    self.logArea.configure(state='disabled')
                    self.logArea.update()

                    genomeName = file2Annotate
                    prot2map = []
                    suffixName = (genomeName.split("/"))[-1]

                    gffFile = open(suffixName+"_annotation.gff","w") #To change from the command line
                    warnFile = open(suffixName+"_annotationWarnings.txt","w") #To change from the command line
                    cdsFile = open(suffixName+"_cds.fasta","w") #To change from the command line
                    protFile = open(suffixName+"_proteins.fasta","w") #To change from the command line


                    os.system("rm -f "+installationDirectory+"annotation/proteinDB/._*")

                    for seq_record in SeqIO.parse(genomeName,"fasta"):
                        genomeSeq = str(seq_record.seq)
                        assemblyName = str(seq_record.id)


                    onlyfiles = [f for f in listdir(installationDirectory+"annotation/proteinDB/") if isfile(join(installationDirectory+"annotation/proteinDB/", f))]
                    for f in onlyfiles:
                        if f[-13:] == "_models.fasta":
                            prot2map.append(f)

                    for f in prot2map:
                        numCodonRefines = 0
                        #Find the best model for the gene *********************************************************
                        locus = f.replace("_models.fasta","")
                        protSeqs = SeqIO.to_dict(SeqIO.parse(installationDirectory+"annotation/proteinDB/"+f,"fasta"))
                        #print "Choosing best match for protein",f
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "\nSearching best match of "+f+" for "+file2AnnotateName+"....\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()

                        os.system(installationDirectory+"resources/blastx -query "+genomeName+" -db "+installationDirectory+"annotation/proteinDB/"+f+" -outfmt 6 | sort   -k 12rn,12rn > outputBlast.txt")
                        blastFile = open("outputBlast.txt")
                        bestProt = ((blastFile.readline()).split("\t"))[1]

                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Found "+bestProt+"\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()   
                        #print "Best match in the following protein:",bestProt

                        sequence = str(protSeqs[bestProt].seq)
                        tempFasta = open("tempFasta.fasta","w")
                        tempFasta.write(">"+locus+"\n"+str(protSeqs[bestProt].seq)+"\n")
                        tempFasta.close()

                        #Run Exonerate on the best model *********************************************************
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Retrieving model....\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update() 
                        os.system(installationDirectory+"resources/exonerate --model protein2genome tempFasta.fasta "+genomeName+" --showtargetgff -s 0 -n 1 --forcegtag --minintron 35  >outputExonerate")
                       
                    
                    
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
                                if not cdsSeq[:3]=="ATG" or not (cdsSeq[:3] =="TTG" and locus=="RL6"): #RL6 start with alternative start codon
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
                                            exon[locus][0]=(int(exon[locus][0][0])-a*3,exon[locus][0][1],exon[locus][0][2])
                                            gene[locus] =(int(exon[locus][0][0])-a*3, int(gene[locus][1]),gene[locus][2])
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
                            gffNote = ""
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
                                                        newExonSet[locus].append((int(item[0]),int(item[0]) + a - newmRNALength -1 ,item[2]))
                                                        exon[locus] = newExonSet[locus]
                                                        gene[locus] = (int(gene[locus][0]),int(item[0]) + a - newmRNALength -1 ,gene[locus][2])
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
                                                    if int(exon[locus][a][1])-int(exon[locus][a][0]) + newmRNALength > a+3:
                                                        newExonSet[locus].append((int(exon[locus][a][1]) - a  + newmRNALength, int(exon[locus][a][1]),exon[locus][a][2]))
                                                        #print "Previous exon locus",exon[locus]
                                                        exon[locus] = newExonSet[locus]
                                                        #print "after exon locus",exon[locus]
                                                        #print gene[locus]
                                                        #gene[locus] = (int(exon[locus][a][1]) - a  + newmRNALength, int(gene[locus][1]), gene[locus][2])
                                                        break

                                                    else:
                                                        newmRNALength = int(exon[locus][a][1]) - int(exon[locus][a][0])
                                                        newExonSet[locus].append((int(exon[locus][a][0]),int(exon[locus][a][1]),exon[locus][a][2]))
                                                cdsSeq = cdsSeq[:a+3]
                                                #print exon[locus][0][1]
                                                exon[locus][0] = (exon[locus][0][0], exon[locus][0][1]-6,exon[locus][0][2])
                                                #print exon[locus][0][1]
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
                                        if item == exon[locus][-1]: #If this is the last exon include the stop codon in the coordinates
                                            gffFile.write(assemblyName+"\texonerate\t"+"CDS"+"\t"+str(item[0])+"\t"+str(int(item[1])+3)+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+"\n")
                                        else:
                                            gffFile.write(assemblyName+"\texonerate\t"+"CDS"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+"\n")
                                        numExon += 1
                                    

                                else: # *********************** Negative Strand
                                    gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(int(gene[locus][0])-3)+"\t"+str(int(gene[locus][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+"\n")
                                    gffFile.write(assemblyName+"\texonerate\t"+"mRNA"+"\t"+str(int(gene[locus][0])-3)+"\t"+str(int(gene[locus][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+"\n")
                                    numExon = 1
                                    for item in exon[locus]:
                                        if item == exon[locus][0]:# and len(exon[locus]) == 1: #If this is the first exon include the stop codon in the coordinates
                                            gffFile.write(assemblyName+"\texonerate\t"+"CDS"+"\t"+str(int(item[0])-3)+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+"\n")
                                        else:
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
                            self.logArea.configure(state='normal')
                            self.logArea.insert(tk.END, "Annotation needs refinement....\n")
                            self.logArea.see(tk.END)
                            self.logArea.configure(state='disabled')
                            self.logArea.update() 
                            print "Annotation needs refinement"
                            exResult.close()
                            #  ***************************************************************************
                            #  ************************* Annotation refinement ***************************
                            #  ***************************************************************************
                            

                            os.system(installationDirectory+"resources/exonerate --model protein2genome tempFasta.fasta "+genomeName+" --showtargetgff -s 0 -n 1 --refine full --forcegtag --minintron 35 >outputExonerate")
                            
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
                                    if not cdsSeq[:3]=="ATG" or not (cdsSeq[:3] =="TTG" and locus=="RL6"): #RL6 start with a non canonical start codon
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
                                                exon[locus][0]=(int(exon[locus][0][0])-a*3,exon[locus][0][1],exon[locus][0][2])
                                                gene[locus] =(int(exon[locus][0][0])-a*3, int(gene[locus][1]),gene[locus][2])
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
                                                        if int(exon[locus][a][1])-int(exon[locus][a][0]) + newmRNALength > a+3:
                                                            newExonSet[locus].append((int(exon[locus][a][1]) - a -3 + newmRNALength, int(exon[locus][a][1]),exon[locus][a][2]))
                                                            exon[locus] = newExonSet[locus]
                                                            gene[locus] = (int(exon[locus][a][1]) - a -3 + newmRNALength, int(gene[locus][1]), gene[locus][2])
                                                            break

                                                        else:
                                                            newmRNALength = int(exon[locus][a][1]) - int(exon[locus][a][0])
                                                            newExonSet[locus].append((int(exon[locus][a][0]),int(exon[locus][a][1]),exon[locus][a][2]))
                                                    cdsSeq = cdsSeq[:a+3]
                                                    exon[locus][0] = (exon[locus][0][0], exon[locus][0][1]-6,exon[locus][0][2])

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
                                            if item == exon[locus][-1]: 
                                                gffFile.write(assemblyName+"\texonerate\t"+"CDS"+"\t"+str(item[0])+"\t"+str(int(item[1])+3)+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+"\n")
                                            else:
                                                gffFile.write(assemblyName+"\texonerate\t"+"CDS"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+"\n")
                                            numExon += 1
                                        

                                    else: # *********************** Negative Strand
                                        gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(int(gene[locus][0])-3)+"\t"+str(int(gene[locus][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+"\n")
                                        gffFile.write(assemblyName+"\texonerate\t"+"mRNA"+"\t"+str(int(gene[locus][0])-3)+"\t"+str(int(gene[locus][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+"\n")
                                        numExon = 1
                                        for item in exon[locus]:
                                            if item == exon[locus][0]:# and len(exon[locus]) == 1: 
                                                gffFile.write(assemblyName+"\texonerate\t"+"CDS"+"\t"+str(int(item[0])-3)+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+"\n")
                                            else:
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






        '''This class configures and populates the toplevel window.
           top is the toplevel containing window.'''
        _bgcolor = '#d9d9d9'  # X11 color: 'gray85'
        _fgcolor = '#000000'  # X11 color: 'black'
        _compcolor = '#d9d9d9' # X11 color: 'gray85'
        _ana1color = '#d9d9d9' # X11 color: 'gray85' 
        _ana2color = '#ececec' # Closest X11 color: 'gray92' 
        self.style = ttk.Style()
        if sys.platform == "win32":
            self.style.theme_use('winnative')
        self.style.configure('.',background=_bgcolor)
        self.style.configure('.',foreground=_fgcolor)
        self.style.configure('.',font="TkDefaultFont")
        self.style.map('.',background=
            [('selected', _compcolor), ('active',_ana2color)])

        w=900
        h=400
        ws = root.winfo_screenwidth() # width of the screen
        hs = root.winfo_screenheight() # height of the screen
        # calculate x and y coordinates for the Tk root window
        x = (ws/2) - (w/2)
        y = (hs/2) - (h/2)

        top.geometry('%dx%d+%d+%d' % (w, h, x, y)) 
        top.title("Annotation tool")
        top.configure(highlightcolor="black")

        self.inputFileLabel = tk.Label(top)
        self.inputFileLabel.configure(text="Input file")
        self.inputFileLabel.place(x=20,y=20,width=70,height=20)

        self.inputFileEntry = tk.Entry(top)
        self.inputFileEntry.place(x=20,y=40,width=750,height=30)
        self.inputFileEntry.insert(0,"Please select a genomes list....")

        self.inputFileButton = tk.Button(top,command=openInputFile)
        self.inputFileButton.place(x=780,y=40,width=100,height=30)
        self.inputFileButton.configure(text="Open file")


        self.outputFolderLabel = tk.Label(top)
        self.outputFolderLabel.configure(text="Output folder")
        self.outputFolderLabel.place(x=20,y=100,width=90,height=20)

        self.outputFolderEntry = tk.Entry(top)
        self.outputFolderEntry.place(x=20,y=120,width=750,height=30)
        self.outputFolderEntry.insert(0,"Please select an output folder....")

        self.outputFolderButton = tk.Button(top,command=openOutputFolder)
        self.outputFolderButton.place(x=780,y=120,width=100,height=30)
        self.outputFolderButton.configure(text="Open folder")



        self.runButton = tk.Button(top,command=mainAnnotationAlgorithm)
        self.runButton.place(x=780,y=300,width=100,height=30)
        self.runButton.configure(text="Run")

        self.runButton = tk.Button(top,command=exitProgram)
        self.runButton.place(x=780,y=350,width=100,height=30)
        self.runButton.configure(text="Exit")

        self.logFileLabel = tk.Label(top)
        self.logFileLabel.place(x=20,y=180,height=20,width=100)
        self.logFileLabel.configure(text="Log window")

        self.logFrame = tk.Frame(top)
        self.logFrame.place(x=20, y=200, height=180, width=750)
        self.logFrame.configure(relief='groove')
        self.logFrame.configure(borderwidth="2")
        self.logArea = tk.Text(top,state='disabled')
        self.logArea.place(x=25,y=205,height=170, width=740)
        self.logArea.configure(background="white",borderwidth=5)
        self.logArea.configure(selectbackground="#c4c4c4")




if __name__ == '__main__':
    vp_start_gui()
