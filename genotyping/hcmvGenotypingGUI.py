#!/usr/bin/python
installationFolder = "/home3/scc20x/Software/mySoftware/GRACy/genotyping"


from pyPdf import PdfFileWriter, PdfFileReader
import StringIO
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter


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
import datetime
import random as rd
from Bio import SeqIO
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Bio.Graphics import GenomeDiagram
from reportlab.lib import colors
from Bio.SeqFeature import SeqFeature, FeatureLocation



def vp_start_gui():
    '''Starting point when module is the main routine.'''
    global val, w, root

    

    root = tk.Tk()
    top = Toplevel1 (root)

    


    #top.InputFileButton.bind('<Button-1>',openInputFile)

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
            inputFile = tkFileDialog.askopenfilename(initialdir = "/",title = "Select an input file")
            inputFile = (inputFile.split("/"))[-1]
            self.InputFileEntry.delete(0,tk.END)
            self.InputFileEntry.insert(0,inputFile)

        def changeDBFile():
            dbfile = tkFileDialog.askopenfilename(initialdir = "/",title = "Select a kmer database")
            dbfile = (dbfile.split("/"))[-1]
            self.dbEntry.delete(0,tk.END)
            self.dbEntry.insert(0,dbfile)

        #**************************************************************************
        #***************** Main Genotyping   Algorithm Start **********************
        #**************************************************************************

        def runGenotyping():

            #Initialize variables
            orderedHyperLoci = ["rl5a","rl6","rl12","rl13","ul1","ul9","ul11","ul20","ul73","ul74","ul120","ul139","ul146"]

            numReads = int(self.numReadsEntry.get())
            dbfile = installationFolder+"/kmerDB/"+self.dbEntry.get()
            NumThreads = int(self.numThreadsEntry.get())
            inputFile = self.InputFileEntry.get()

            datasets = []
            infile=open(inputFile)
            while True:
                read1 = infile.readline().rstrip()
                if not read1:
                    break  
                read2 = infile.readline().rstrip()
                if not read2:
                    break

                datasets.append((read1,read2))
            
            numDatasets = len(datasets)
            progressBarIncrement = 100/(numDatasets*20)
            step = 0
            for dataset in datasets:
                read1 = dataset[0]
                read2 = dataset[1]

                #Extract path and files info
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
                logFile = open(fileRoot1+"_logFile.txt","w")
                now = datetime.datetime.now()
                logFile.write("Genotyping sample "+fileRoot1+"  started at "+now.strftime("%H:%M")+"\n")
                
                fileList = "list"+str(rd.randint(0,1000000))+".txt"
                listFile = open(fileList,"w")
                listFile.write(read1+"\n"+read2+"\n")
                listFile.close()

                #Perform deduplication
                now = datetime.datetime.now()
                logFile.write("Deduplication for dataset "+fileRoot1+"  started at "+now.strftime("%H:%M")+"\n")

                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "******************************************************************************\n")
                self.logArea.insert(tk.END, "Genotyping dataset "+fileRoot1+"/"+fileRoot2+"\n")
                self.logArea.insert(tk.END, "******************************************************************************\n\n")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()


                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "Performing deduplication....")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()

                #print "Performing deduplication...."
                dedupFile1 = fileRoot1+"_noDup_1.fq"
                dedupFile2 = fileRoot2+"_noDup_2.fq"
                os.system("fastuniq -i "+fileList+" -t q -o "+dedupFile1+" -p "+dedupFile2)
                os.system("rm -f "+fileList)

                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "Done!\n")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()

                now = datetime.datetime.now()
                logFile.write("Deduplication for dataset "+fileRoot1+"  ended at "+now.strftime("%H:%M")+"\n\n")
                step = step+1
                self.progressbar['value']=int( step*progressBarIncrement)
                self.progressbar.update()

                #******************* Calculate average coverage for deduplicated reads   **************
                now = datetime.datetime.now()
                logFile.write("Reference coverage calculation for dataset "+fileRoot1+"  started at "+now.strftime("%H:%M")+"\n")

                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "Performing reference alignment....")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()
                os.system("bowtie2 -1 "+dedupFile1+" -2 "+dedupFile2+" -x "+installationFolder+"/fastaFiles/hcmvReference -S alignmenthsbfy43223.sam >null 2>&1")
                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "Done!\n")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()
                step = step+1
                self.progressbar['value']=int( step*progressBarIncrement)
                self.progressbar.update()
                
                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "Converting sam to bam....")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()
                os.system("samtools view -bS -h alignmenthsbfy43223.sam > alignmenthsbfy43223.bam")
                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "Done!\n")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()
                step = step+1
                self.progressbar['value']=int( step*progressBarIncrement)
                self.progressbar.update()
                
                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "Sorting bam....")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()
                os.system("samtools sort -o alignmenthsbfy43223_sorted.bam alignmenthsbfy43223.bam")
                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "Done!\n")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()
                step = step+1
                self.progressbar['value']=int( step*progressBarIncrement)
                self.progressbar.update()

                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "Calculating average coverage....")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()
                os.system("samtools depth  alignmenthsbfy43223_sorted.bam  |  awk '{sum+=$3} END { print sum/NR}' >avCoverage.txt")
                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "Done!\n")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()
                step = step+1
                self.progressbar['value']=int( step*progressBarIncrement)
                self.progressbar.update()

                avCovFile = open("avCoverage.txt")
                avCov = float(avCovFile.readline().rstrip())
                avCovFile.close()
                detectionTreshold = float(avCov*float(self.CutoffText.get()))  
                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "Average coverage for deduplicated reads: "+str(avCov)+"\n")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()
                os.system("rm -f alignmenthsbfy43223.* avCoverage.txt "+fileList)
                now = datetime.datetime.now()
                logFile.write("Reference coverage calculation for dataset "+fileRoot1+"  ended at "+now.strftime("%H:%M")+"\n\n")

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
                                kmerLength = len(item)
                                geneKmers[(fields[0],fields[1])].append(item)


                #Get sequences in memory
                reads = []
                numSeq = 0
                overallGeneInfo = {}
                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "Caclculate kmer frequencies for file "+dedupFile1+"....")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()
                os.system("jellyfish count -m 17 -s 100M -t 8 -C "+dedupFile1+" -o "+dedupFile1+"_kmerCount.jf")

                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "Done!\n")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()
                step = step+1
                self.progressbar['value']=int( step*progressBarIncrement)
                self.progressbar.update()

                numSeq = 0
                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "Caclculate kmer frequencies for file "+dedupFile2+"....")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()

                os.system("jellyfish count -m 17 -s 100M -t 8 -C "+dedupFile2+" -o "+dedupFile2+"_kmerCount.jf")

                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "Done!\n")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()
                

                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "Merging kmer files....")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()

                os.system("jellyfish merge "+dedupFile1+"_kmerCount.jf "+dedupFile2+"_kmerCount.jf")
                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "Done!\n")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()

                step = step+1
                self.progressbar['value']=int( step*progressBarIncrement)
                self.progressbar.update()

                readsKmer = {}

                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "Loading kmers for sequences in "+dedupFile1+" into memory....")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()
                for seq_record in SeqIO.parse(dedupFile1,"fastq"):
                    sequence = str(seq_record.seq)
                    revSequence = str(seq_record.seq.reverse_complement())
                    reads.append(str(seq_record.seq))
                    reads.append(str(seq_record.seq.reverse_complement()))
                    for a in range(0,len(sequence)-kmerLength+1):
                        if not sequence[a:a+kmerLength] in readsKmer:
                            readsKmer[sequence[a:a+kmerLength]] = []
                        readsKmer[sequence[a:a+kmerLength]].append(str(seq_record.id))

                    for a in range(0,len(revSequence)-kmerLength+1):
                        if not revSequence[a:a+kmerLength] in readsKmer:
                            readsKmer[revSequence[a:a+kmerLength]] = []
                        readsKmer[revSequence[a:a+kmerLength]].append(str(seq_record.id))

                    numSeq +=1
                    if numSeq == int(numReads/2):
                        break

                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "Done!\n")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()


                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "Loading kmers for sequences in "+dedupFile2+" into memory....")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()

                for seq_record in SeqIO.parse(dedupFile2,"fastq"):
                    sequence = str(seq_record.seq)
                    revSequence = str(seq_record.seq.reverse_complement())
                    reads.append(str(seq_record.seq))
                    reads.append(str(seq_record.seq.reverse_complement()))
                    for a in range(0,len(sequence)-16):
                        if not sequence[a:a+17] in readsKmer:
                            readsKmer[sequence[a:a+17]] = []
                        readsKmer[sequence[a:a+17]].append(str(seq_record.id))

                    for a in range(0,len(revSequence)-16):
                        if not revSequence[a:a+17] in readsKmer:
                            readsKmer[revSequence[a:a+17]] = []
                        readsKmer[revSequence[a:a+17]].append(str(seq_record.id))

                    numSeq +=1
                    if numSeq == int(numReads/2):
                        break

                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "Done!\n")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()


                now = datetime.datetime.now()
                logFile.write("Genotyping for dataset "+fileRoot1+"  started at "+now.strftime("%H:%M")+"\n")

                #Start searching for signatures
                for gene in orderedHyperLoci:
                    self.logArea.configure(state='normal')
                    self.logArea.insert(tk.END, "Genotyping gene "+gene+"....")
                    self.logArea.see(tk.END)
                    self.logArea.configure(state='disabled')
                    self.logArea.update()
                    
                    matchedReads = {}
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
                    numMatchedKmers = {}
                    for gr in specificKmerGroup:


                        if not gr in countSeq:
                            countSeq[gr] = 0

                        if not gr in matchedReads:
                            matchedReads[gr] = set()

                        if not gr in numMatchedKmers:
                            numMatchedKmers[gr] = 0
                        
                        command = "jellyfish query mer_counts_merged.jf "
                        for querySeq in specificKmerGroup[gr]:
                            command += querySeq
                            command += " "
                        command += " >counts.txt"

                        os.system(command)

                        countFile = open("counts.txt")
                        
                        
                        while True:
                            countLine = countFile.readline().rstrip()
                            if not countLine:
                                break
                            countFields = countLine.split(" ")
                            if int(countFields[1])>0:
                                numMatchedKmers[gr] +=1
                                for item in readsKmer[countFields[0]]:
                                    matchedReads[gr].add(item)
                        countFile.close()
                        logFile.write(gene+"\t"+"For genotype "+gr+" there are "+str(len(matchedReads[gr]))+" reads matching\n")
                        if len(matchedReads[gr])>=detectionTreshold:
                            countSeq[gr] = len(matchedReads[gr])
                            totCount = totCount + len(matchedReads[gr])

                    self.logArea.configure(state='normal')
                    self.logArea.insert(tk.END, "Done!\n")
                    self.logArea.see(tk.END)
                    self.logArea.configure(state='disabled')
                    self.logArea.update()
                    
                    outfile.write(gene)
                    
                    
                    for gr in countSeq:
                        if countSeq[gr]>0:
                            percentage = float(countSeq[gr])/float(totCount)
                            outfile.write("\t"+gr+"\t"+str(percentage)+"\t"+str(totCount)+"\t"+str(numMatchedKmers[gr]))
                            self.logArea.configure(state='normal')
                            self.logArea.insert(tk.END, "Found signature for genotype "+gr+" with a percentage of "+(str(percentage*100))[:5]+"%\n")
                            self.logArea.see(tk.END)
                            self.logArea.configure(state='disabled')
                            self.logArea.update()

                    outfile.write("\n")
                    step = step+1
                    self.progressbar['value']=int( step*progressBarIncrement)
                    self.progressbar.update()

                now = datetime.datetime.now()
                logFile.write("Genotyping sample "+fileRoot1+"  ended at "+now.strftime("%H:%M")+"\n")
                logFile.close()
                outfile.close()
                os.system("rm -f "+dedupFile1+" "+dedupFile2+" *.jf")

        #**************************************************************************
        #***************** Main Genotyping   Algorithm End * **********************
        #**************************************************************************



        #**************************************************************************
        #***************** Main Plotting Algorithm Start  *************************
        #**************************************************************************
        def plotGenotypes():
            packet = StringIO.StringIO()
            inputFiles = tkFileDialog.askopenfilenames(initialdir = "./",title = "Select an input file")
            
            
            #Initialize variables
            orderedHyperLoci = ["rl5a","rl6","rl12","rl13","ul1","ul9","ul11","ul20","ul73","ul74","ul120","ul139","ul146"]
            genotypeColors = {}
            allGenotypes = ["G1","G10","G11","G12","G13","G14","G1A","G1B","G1C","G2","G2A","G2B","G3","G3A","G3B","G4","G4A","G4B","G4C","G4D","G5","G6","G7","G8","G9"]
            
            legColour = ["#8b8989","#fffaf0","#faebd7","#eed5b7","#9932cc","#cdcdc1","#c1cdc1","#32cd32","#ffe4e1","#2f4f4f","#778899","#000080","#6495ed","#8470ff","#00bfff","#00ffff","#66cdaa","#006400","#8fbc8f","#7fff00","#cd5c5c","#a0522d","#f4a460","#b22222","#ff0000"]
            #Assign colour to genotypes
            rd.seed(1200)
            colorStep = 0
            for g in allGenotypes:
                colorStep +=1
                if not g in genotypeColors:
                    r = lambda: rd.randint(0,100)
                    genotypeColors[g] = legColour[colorStep-1]
                    #genotypeColors[g]=(float(r())/100.0,float(r())/100.0,float(r())/100)
                    #genotypeColors[g]=(float(colorStep)/100.0,float(90-colorStep+7)/100.0,float(colorStep+21)/100)
            
            #Plot legend
            numLeg = 0
            for gr in genotypeColors:
                numLeg += 1
                plt.rc('figure', figsize=(5, 5))
                
                plt.subplot(len(genotypeColors),1,numLeg)
                plt.axis('equal')
                plt.pie([100],labels=[gr],colors=[genotypeColors[gr]])
                plt.axis('equal')
                #plt.text(-5,0,"Ciao2",verticalalignment='center',horizontalalignment='left')
                if numLeg == 1:
                    plt.title("Colour legend")

            plt.savefig('legend.pdf',format='pdf',dpi=600, bbox_inches = 'tight')



            #Print hypervariable genes ideogram
            numHyper = len(orderedHyperLoci)
            maxLen = numHyper*1200
            gd_diagram = GenomeDiagram.Diagram("hcmv genotyping")
            gd_track_for_features = gd_diagram.new_track(1, name="Ideogram",height=0.02,scale=0)
            gd_feature_set = gd_track_for_features.new_set()
            for a in range(len(orderedHyperLoci)):
                feature = SeqFeature(FeatureLocation(1200*a+100, 1200*a+1000+100), type=orderedHyperLoci[a],strand=+1)
                gd_feature_set.add_feature(feature, color=colors.gray, label=True,label_size=15, label_angle=-45,label_position="middle")


            inputDirectory = "/".join(inputFiles[0].split("/")[:-1])

            #Collect data to plot
            samplesToPlot = {}
            for sample in inputFiles:
                sampleName = (sample.split("/"))[-1]
                if not sampleName in samplesToPlot:
                    samplesToPlot[sampleName] = {}
            sampleList = []
            for sample in inputFiles:
                sampleList.append(sample.split("/")[-1]) 
                sampleName = (sample.split("/"))[-1]
                genotypes = open(sample)
                while True:
                    line = genotypes.readline().rstrip()
                    if not line:
                        break
                    fields = line.split("\t")
                    if not fields[0] in samplesToPlot[sampleName]:
                        samplesToPlot[sampleName][fields[0]] = [] 
                        for a in range(1,len(fields)-1,+4):
                            samplesToPlot[sampleName][fields[0]].append(fields[a])
                            samplesToPlot[sampleName][fields[0]].append(fields[a+1])
                            print fields[a],fields[a+1]

                        #for item in fields[1:]:
                        #    samplesToPlot[sample][fields[0]].append(item)
                genotypes.close()
                

            
            print "samplesToPlot"
            print samplesToPlot
            print "sampleList"
            print sampleList
            #Add plot to final diagram
            for sample in sampleList:
                gd_track_for_genotypes = gd_diagram.new_track(1, name=sample,height=0.2, scale=0)
                gd_genotype_set = gd_track_for_genotypes.new_set()
                for a in range(len(orderedHyperLoci)):
                    start = a*1200
                    for b in range(0,len(samplesToPlot[sample][orderedHyperLoci[a]]),+2):
                        #print orderedHyperLoci[a], samplesToPlot[sample][orderedHyperLoci[a]]

                        feature = SeqFeature(FeatureLocation(start+100, start+100+int(float(samplesToPlot[sample][orderedHyperLoci[a]][b+1])*1000) ), type=samplesToPlot[sample][orderedHyperLoci[a]][b],strand=+1)
                        gd_genotype_set.add_feature(feature,color = genotypeColors[    samplesToPlot[sample][orderedHyperLoci[a]][b] ])
                        #print start+100, start+100+int(float(samplesToPlot[sample][orderedHyperLoci[a]][b+1])*1000)
                        start=start+int(float(samplesToPlot[sample][orderedHyperLoci[a]][b+1])*1000)




            #Print plot on file
            gd_track_for_fake = gd_diagram.new_track(1, name="first genotype",height=0.6,scale=0)
            gd_diagram.draw(format="cicular",orientation="landscape", pagesize='A4', start=0, end=maxLen, circle_core=0.4, tracklines=0, track_size=1,x=0.2)
            gd_diagram.write("plot.pdf", "PDF",dpi=300)

            
            can = canvas.Canvas(packet, pagesize=letter)
            startY = 400
            for sample in sampleList:
                startY -= 20
                can.drawString(300,startY, sample)
            can.save()
            packet.seek(0)
            new_pdf = PdfFileReader(packet)
            existing_pdf = PdfFileReader(file("plot.pdf", "rb"))
 
            output = PdfFileWriter()
            page = existing_pdf.getPage(0)
            page.mergePage(new_pdf.getPage(0))
            output.addPage(page)
            outputStream = file("plot2.pdf", "wb")
            output.write(outputStream)
            outputStream.close()

            can = None
            new_pdf = None
            existing_pdf = None
            output = None
            page = None
            outputStream = None

            





        def exitProgram():
            exit()

        
            

        


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

        top.geometry("960x460+609+943")
        top.title("Genotyping")
        top.configure(highlightcolor="black")


        self.InputFileLabel = tk.Label(top)
        self.InputFileLabel.place(x=20, y=20,height=20, width=60)
        self.InputFileLabel.configure(text="Input file")

        self.InputFileEntry = tk.Entry(top)
        self.InputFileEntry.place(x=20, y=40,height=30, width=300)
        self.InputFileEntry.insert(0,"/home3/scc20x/Software/mySoftware/GRACy/genotyping/inputFile.txt")

        self.InputFileButton = tk.Button(top,command=openInputFile)
        self.InputFileButton.place(x=340,y=40,height=30, width=100)
        self.InputFileButton.configure(text="Open file")

        self.dblabel = tk.Label(top)
        self.dblabel.place(x=480, y=20,height=20, width=100)
        self.dblabel.configure(text="kmer database")

        self.dbEntry = tk.Entry(top)
        self.dbEntry.place(x=480, y=40,height=30, width=200)
        self.dbEntry.insert(0,"mainDB_seqs_filtered.txt")

        self.changeDB = tk.Button(top)
        self.changeDB.place(x=700, y=40,height=30, width=100)
        self.changeDB.configure(text="Change DB")

        self.newDB = tk.Button(top)
        self.newDB.place(x=820, y=40,height=30, width=100)
        self.newDB.configure(text="New DB")
        
        self.numReadsLabel = tk.Label(top)
        self.numReadsLabel.place(x=20, y=100, height=20, width=120)
        self.numReadsLabel.configure(text="Number of reads")

        self.numReadsEntry = tk.Entry(top, justify='right')
        self.numReadsEntry.place(x=20, y=120,height=30, width=150)
        self.numReadsEntry.insert(0,"10000000000000")

        self.CutoffLabel = tk.Label(top)
        self.CutoffLabel.place(x=200, y=100,height=20, width=120)
        self.CutoffLabel.configure(text="Detection treshold")

        self.CutoffText = tk.Entry(top, justify='right')
        self.CutoffText.place(x=200, y=120,height=30, width=150)
        self.CutoffText.insert(0,"0.2")

        self.logFileLabel = tk.Label(top)
        self.logFileLabel.place(x=20,y=200,height=20,width=100)
        self.logFileLabel.configure(text="Log window")

        self.logFrame = tk.Frame(top)
        self.logFrame.place(x=20, y=220, height=180, width=920)
        self.logFrame.configure(relief='groove')
        self.logFrame.configure(borderwidth="2")
        self.logFrame.configure(relief='groove')
        self.logFrame.configure(width=125)
        self.logArea = tk.Text(top,state='disabled')
        self.logArea.place(x=25,y=225,height=170, width=910)
        self.logArea.configure(background="white",borderwidth=5)
        self.logArea.configure(selectbackground="#c4c4c4")

        self.progressbar=ttk.Progressbar(top,orient="horizontal",length=920,mode="determinate")
        self.progressbar.place(x=20,y=420)
        self.progressbar['maximum'] = 100

        self.numThreadsLabel = tk.Label(top)
        self.numThreadsLabel.place(x=380,y=100,height=20, width=100)
        self.numThreadsLabel.configure(text="Num. threads")

        self.numThreadsEntry = tk.Entry(top, justify='right')
        self.numThreadsEntry.place(x=380, y=120,height=30, width=150)
        self.numThreadsEntry.insert(0,"8")

        self.runButton = tk.Button(top,command=runGenotyping)
        self.runButton.place(x=560,y=120,height=30,width=130)
        self.runButton.configure(text="Genotype!")

        self.plotButton = tk.Button(top,command=plotGenotypes)
        self.plotButton.place(x=700,y=120,height=30,width=130)
        self.plotButton.configure(text="Plot")

        self.exitButton = tk.Button(top,command=exitProgram)
        self.exitButton.place(x=840,y=120,height=30,width=80)
        self.exitButton.configure(text="Exit")






if __name__ == '__main__':
    vp_start_gui()