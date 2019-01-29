#!/usr/bin/python
installationFolder = "/home3/scc20x/Software/mySoftware/GRACy/genotyping"


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
        #***************** Main Algorithm Start ***********************************
        #**************************************************************************

        def runGenotyping():

            #Initialize variables
            orderedHyperLoci = ["rl5a","rl6","rl12","rl13","ul1","ul9","ul11","ul20","ul73","ul74","ul120","ul139","ul146"]
            allGroups = set()
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
            for dataset in datasets:
                step = 0
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

                print "Performing deduplication...."
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
                detectionTreshold = int(avCov*float(self.CutoffText.get()))  
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
                                geneKmers[(fields[0],fields[1])].append(item)


                #Get sequences in memory
                reads = []
                numSeq = 0
                overallGeneInfo = {}
                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "Loading reads from "+fileRoot1+"  into memory....")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()
                for seq_record in SeqIO.parse(read1,"fastq"):
                    reads.append(str(seq_record.seq))
                    reads.append(str(seq_record.seq.reverse_complement()))
                    numSeq +=1
                    if numSeq == int(numReads/2):
                        break
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
                self.logArea.insert(tk.END, "Loading reads from "+fileRoot2+"  into memory....")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()
                for seq_record in SeqIO.parse(read2,"fastq"):
                    reads.append(str(seq_record.seq))
                    reads.append(str(seq_record.seq.reverse_complement()))
                    numSeq +=1
                    if numSeq == int(numReads/2):
                        break
                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "Done!\n")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()
                step = step+1
                self.progressbar['value']=int( step*progressBarIncrement)
                self.progressbar.update()



                now = datetime.datetime.now()
                logFile.write("Genotyping for dataset "+fileRoot1+"  started at "+now.strftime("%H:%M")+"\n")

                #Start searching for signatures
                for gene in orderedHyperLoci:
                    self.logArea.configure(state='normal')
                    self.logArea.insert(tk.END, "Genotyping gene "+gene+"....")
                    self.logArea.see(tk.END)
                    self.logArea.configure(state='disabled')
                    self.logArea.update()
                    
                    
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
                        if countSeq[gr]< detectionTreshold:
                            totCount = totCount - countSeq[gr]
                            countSeq[gr]=0
                        if len(matchingKmer)>0:
                            logFile.write("For group "+gr+" there are "+str(len(matchingKmer))+" matched kmers over a total of "+str(len(geneKmers[(gene,gr)]))+". Number of matches: "+str(countSeq[gr]))
                                 
                    self.logArea.configure(state='normal')
                    self.logArea.insert(tk.END, "Done!\n")
                    self.logArea.see(tk.END)
                    self.logArea.configure(state='disabled')
                    self.logArea.update()
                    averageCoverage = {}
                    outfile.write(gene)
                    for gr in countSeq:
                        if not gr in averageCoverage and totCount>0:
                            #print gene,gr,float(countSeq[gr]),float(totCount),detectionTreshold
                            averageCoverage[gr] = float(countSeq[gr])/float(totCount)
                            if averageCoverage[gr] >0.02: #Here as treshold may be inserted
                                #print "Values",float(countSeq[gr]),float(totCount)
                                outfile.write("\t"+gr+"\t"+str(averageCoverage[gr]))
                                self.logArea.configure(state='normal')
                                self.logArea.insert(tk.END, "Found signature for genotype "+gr+" with a percentage of "+str(averageCoverage[gr]*100)+"%\n")
                                self.logArea.see(tk.END)
                                self.logArea.configure(state='disabled')
                                self.logArea.update()

                            #print "Statistics2",gr,averageCoverage2[gr],countSeq[gr],totCount

                    outfile.write("\n")
                    step = step+1
                    self.progressbar['value']=int( step*progressBarIncrement)
                    self.progressbar.update()

                now = datetime.datetime.now()
                logFile.write("Genotyping sample "+fileRoot1+"  ended at "+now.strftime("%H:%M")+"\n")
                logFile.close()
                outfile.close()

        #**************************************************************************
        #***************** Main Algorithm End *************************************
        #**************************************************************************

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
        self.dbEntry.insert(0,"mainDB_seqs.txt")

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
        self.numReadsEntry.place(x=20, y=120,height=30, width=200)
        self.numReadsEntry.insert(0,"1000000")

        self.CutoffLabel = tk.Label(top)
        self.CutoffLabel.place(x=260, y=100,height=20, width=120)
        self.CutoffLabel.configure(text="Detection treshold")

        self.CutoffText = tk.Entry(top, justify='right')
        self.CutoffText.place(x=260, y=120,height=30, width=200)
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
        self.numThreadsLabel.place(x=480,y=100,height=20, width=100)
        self.numThreadsLabel.configure(text="Num. threads")

        self.numThreadsEntry = tk.Entry(top, justify='right')
        self.numThreadsEntry.place(x=480, y=120,height=30, width=200)
        self.numThreadsEntry.insert(0,"8")

        self.runButton = tk.Button(top,command=runGenotyping)
        self.runButton.place(x=700,y=120,height=30,width=100)
        self.runButton.configure(text="Run")

        self.exitButton = tk.Button(top,command=exitProgram)
        self.exitButton.place(x=820,y=120,height=30,width=100)
        self.exitButton.configure(text="Exit")






if __name__ == '__main__':
    vp_start_gui()