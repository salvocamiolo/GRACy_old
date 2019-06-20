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


import sys,os
import Tkinter, Tkconstants, tkFileDialog

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
        
        
        def openReferenceFile():
            referenceFile = tkFileDialog.askopenfilename(initialdir = "./",title = "Select reference file")
            self.referenceFileEntry.delete(0,tk.END)
            self.referenceFileEntry.insert(0,referenceFile)
            self.referenceFileEntry.icursor(100)
        
        def openCDSFile():
            cdsFile = tkFileDialog.askopenfilename(initialdir = "./",title = "Select CDS file")
            self.cdsFileEntry.delete(0,tk.END)
            self.cdsFileEntry.insert(0,cdsFile)

        def openGFFFile():
            gffFile = tkFileDialog.askopenfilename(initialdir = "./",title = "Select GFF file")
            self.gffFileEntry.delete(0,tk.END)
            self.gffFileEntry.insert(0,gffFile)

        def openPos2PlotFile():
            pos2plotFile = tkFileDialog.askopenfilename(initialdir = "./",title = "Select file")
            self.posToPlotEntry.delete(0,tk.END)
            self.posToPlotEntry.insert(0,pos2plotFile)
        
        def openInputFile():
            inputFile = tkFileDialog.askopenfilename(initialdir = "./",title = "Select input file")
            self.inputFileEntry.delete(0,tk.END)
            self.inputFileEntry.insert(0,inputFile)

        def openOutputFolder():
            outputolder = tkFileDialog.askdirectory(initialdir = "./",title = "Select folder")
            self.outputFolderEntry.delete(0,tk.END)
            self.outputFolderEntry.insert(0,outputolder)


        def mainAlgorithm():
            if self.referenceFileEntry.get() == "Please select a file...." or self.inputFileEntry.get()=="Please select a file...." or self.gffFileEntry.get() =="Please select a file...." or self.cdsFileEntry.get()=="Please select a file...." or self.outputFolderEntry.get()=="Please select a folder....":
                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "Please enter a valid input, reference, cds, annotation files and output folder  to continue!\n")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()

            else:
                infile = open(self.inputFileEntry.get())
                line = infile.readline().rstrip()
                suffix = (line.split("."))[-1]
                numThreads = self.numThreadsEntry.get()
                outputFolder = self.outputFolderEntry.get()
                referenceFile = self.referenceFileEntry.get()
                if suffix == "fastq" or suffix == "fq":
                    self.logArea.configure(state='normal')
                    self.logArea.insert(tk.END, "The input file contains fastq reads!\n")
                    self.logArea.see(tk.END)
                    self.logArea.configure(state='disabled')
                    self.logArea.update()


                    #***********************************************************************
                    #********************* Main algorithm reads alignment start ************
                    #***********************************************************************

                    infile.close()
                    
                    os.system(installationDirectory+"resources/bowtie2-build "+referenceFile+" reference -q")
                    infile = open(self.inputFileEntry.get())

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
                        
                        print "Analising sample",sampleName
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Analising sample"+sampleName+"....\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()

                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "*  Reads quality filtering before alignment\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()

                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "*  Running Trimgalore....\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()

                        

                        print "Reads quality filtering before alignment"
                        print "Running Trimgalore"
                        prefix1 = ((read1.split("/"))[-1]).replace(".fastq","")
                        prefix2 = ((read2.split("/"))[-1]).replace(".fastq","")
                        os.system(installationDirectory+"resources/trim_galore --paired -q 30  "+read1+" "+read2)

                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "*  Performing deduplication....\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()


                        print "Performing deduplication"
                        os.system("echo "+prefix1+"_val_1.fq > inputFastUniq")
                        os.system("echo "+prefix2+"_val_2.fq >> inputFastUniq")
                        os.system(installationDirectory+"resources/fastuniq -i inputFastUniq -t q -o trimmed_dedup_1.fastq -p trimmed_dedup_2.fastq")

                        print "Performing prinseq quality filtering"

                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "*  Running Prinseq....\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        os.system(installationDirectory+"resources/prinseq -fastq trimmed_dedup_1.fastq  -fastq2 trimmed_dedup_2.fastq -min_qual_mean 25 -trim_qual_right 30  -trim_qual_window 15 -trim_qual_step 5 -min_len 80 -out_bad null -out_good trimmed_dedup_pr")


                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "*  Reads alignment on reference\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()

                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "*  *  Performing alignment....\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()

                        print "Aigning reads to reference"
                        os.system(installationDirectory+"resources/bowtie2 --end-to-end -1 trimmed_dedup_pr_1.fastq -2 trimmed_dedup_pr_2.fastq -x reference -S alignment.sam -p "+numThreads)
                        
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "*  *  Converting sam to bam....\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        
                        print "Converting sam to bam"
                        os.system(installationDirectory+"resources/samtools view -bS -h alignment.sam >alignment.bam")

                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "*  *  Sorting bam....\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()

                        print "Sorting bam"
                        os.system(installationDirectory+"resources/samtools sort -o "+sampleName+"_alignment_sorted.bam alignment.bam")


                        print "Calling snps with lofreq"
                        #Analyze snps with lowfreq
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "*  Calling SNPs with lofreq\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        os.system(installationDirectory+"resources/lofreq call -f "+referenceFile+" -o "+sampleName+"_SNPs.vcf "+sampleName+"_alignment_sorted.bam")
                        
                        #Analyze indes using the GATK pipeline
                        
                        os.system(installationDirectory+"resources/lofreq call -f "+referenceFile+" -o "+sampleName+".vcf "+sampleName+"_alignment_sorted.bam")
                        
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "*  Calling INDELs with GATK\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "*  *  Extracting mapped reads\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        os.system(installationDirectory+"resources/samtools view -bF 4 "+sampleName+"_alignment_sorted.bam >mapped.bam")
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "*  *  Adding group names....\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()

                        os.system("java -jar -XX:ParallelGCThreads="+numThreads+" "+installationDirectory+"resources/picard.jar AddOrReplaceReadGroups I=mapped.bam O=rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=Ilumina RGPU=machine RGSM=Consensus")
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "*  *  Deduplicating....\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        os.system("java -jar -XX:ParallelGCThreads"+numThreads+" "+installationDirectory+"resources/picard.jar MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics")
                        os.system("java -jar "+installationDirectory+"resources/picard.jar CreateSequenceDictionary R="+referenceFile)
                        
                        os.system(installationDirectory+"resources/samtools faidx "+referenceFile)

                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "*  *  Calling INDELs....\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        os.system("java -jar  "+installationDirectory+"resources/GenomeAnalysisTK.jar -T  HaplotypeCaller -R "+referenceFile+" -I dedupped.bam  -o output.vcf -A StrandAlleleCountsBySample")
                        os.system("mv output.vcf "+sampleName+"_indels.vcf")
                            
                        
                        os.system("mv "+sampleName+"_alignment_sorted.bam "+outputFolder+"/")
                        os.system("mv *.vcf "+outputFolder+"/")
                        os.system("rm alignment* -f")
                        os.system("rm -f inputFastUniq *_val_1.fq *_val_2.fq trimmed_dedup_*")
                        os.system(" ls "+outputFolder+"/*SNPs.vcf >vcfFilesToExamine")
                        os.system("python "+installationDirectory+"snpAnalysis/polyAn.py vcfFilesToExamine"+" "+self.cdsFileEntry.get()+" "+self.gffFileEntry.get())
                        os.system("rm vcfFilesToExamine *.dict *.fai mapped.bam *.bt2 rg_added_sorted.bam *trimming_report.txt")
                        os.system("ls *_snpEffect.txt > file2Plot")
                        os.system("python "+installationDirectory+"snpAnalysis/buildTable.py file2Plot")
                        os.system("mv  *_snpEffect.txt *_snpFreq.txt "+outputFolder+"/" )
                        

                if suffix == "vcf":
                    self.logArea.configure(state='normal')
                    self.logArea.insert(tk.END, "The input file contains vcf datasets!\n")
                    self.logArea.see(tk.END)
                    self.logArea.configure(state='disabled')
                    self.logArea.update()
                    os.system("python "+installationDirectory+"snpAnalysis/polyAn.py "+self.inputFileEntry.get()+" "+self.cdsFileEntry.get()+" "+self.gffFileEntry.get())
                    os.system("rm vcfFilesToExamine *.dict *.fai mapped.bam *.bt2 rg_added_sorted.bam *trimming_report.txt")
                    os.system("ls *_snpEffect.txt > file2Plot")
                    os.system("python "+installationDirectory+"snpAnalysis/buildTable.py file2Plot")
                    os.system("mv  *_snpEffect.txt *_snpFreq.txt "+outputFolder+"/" )
                    

                    

                if not suffix == "fastq" and not suffix == "vcf" and not suffix == "fq":
                    self.logArea.configure(state='normal')
                    self.logArea.insert(tk.END, "The input file does not contain datasets GRACy can recognise!\n")
                    self.logArea.see(tk.END)
                    self.logArea.configure(state='disabled')
                    self.logArea.update()



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

        w=900
        h=510
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
        self.inputFileEntry.insert(0,"Please select a file....")

        self.inputFileButton = tk.Button(top,command=openInputFile)
        self.inputFileButton.place(x=780,y=40,width=100,height=30)
        self.inputFileButton.configure(text="Open file")

        #self.posToPlotLabel = tk.Label(top)
        #self.posToPlotLabel.configure(text="Positions to plot")
        #self.posToPlotLabel.place(x=470,y=20,width=110,height=20)

        #self.posToPlotEntry = tk.Entry(top)
        #self.posToPlotEntry.place(x=470,y=40,width=300,height=30)
        #self.posToPlotEntry.insert(0,"Please select a file....")

        #self.posToPlotButton = tk.Button(top,command=openPos2PlotFile)
        #self.posToPlotButton.place(x=780,y=40,width=100,height=30)
        #self.posToPlotButton.configure(text="Open file")


        self.cdsFileLabel = tk.Label(top)
        self.cdsFileLabel.configure(text="CDS file")
        self.cdsFileLabel.place(x=20,y=80,width=60,height=20)

        self.cdsFileEntry = tk.Entry(top)
        self.cdsFileEntry.place(x=20,y=100,width=300,height=30)
        self.cdsFileEntry.insert(0,"Please select a file....")

        self.cdsFileButton = tk.Button(top,command=openCDSFile)
        self.cdsFileButton.place(x=330,y=100,width=100,height=30)
        self.cdsFileButton.configure(text="Open file")


        self.gffFileLabel = tk.Label(top)
        self.gffFileLabel.configure(text="Annotation file")
        self.gffFileLabel.place(x=470,y=80,width=110,height=20)

        self.gffFileEntry = tk.Entry(top)
        self.gffFileEntry.place(x=470,y=100,width=300,height=30)
        self.gffFileEntry.insert(0,"Please select a file....")

        self.gffFileButton = tk.Button(top,command=openGFFFile)
        self.gffFileButton.place(x=780,y=100,width=100,height=30)
        self.gffFileButton.configure(text="Open file")


        self.referenceFileLabel = tk.Label(top)
        self.referenceFileLabel.configure(text="Reference file")
        self.referenceFileLabel.place(x=20,y=140,width=105,height=20)

        self.referenceFileEntry = tk.Entry(top)
        self.referenceFileEntry.place(x=20,y=160,width=300,height=30)
        self.referenceFileEntry.insert(0,"Please select a file....")

        self.referenceFileButton = tk.Button(top,command=openReferenceFile)
        self.referenceFileButton.place(x=330,y=160,width=100,height=30)
        self.referenceFileButton.configure(text="Open file")


        self.outputFolderLabel = tk.Label(top)
        self.outputFolderLabel.configure(text="Output folder")
        self.outputFolderLabel.place(x=470,y=140,width=105,height=20)

        self.outputFolderEntry = tk.Entry(top)
        self.outputFolderEntry.place(x=470,y=160,width=300,height=30)
        self.outputFolderEntry.insert(0,"Please select a folder....")

        self.outputFolderButton = tk.Button(top,command=openOutputFolder)
        self.outputFolderButton.place(x=780,y=160,width=100,height=30)
        self.outputFolderButton.configure(text="Open folder")


        #self.freqCutoffLabel = tk.Label(top)
        #self.freqCutoffLabel.configure(text="Freq cutoff")
        #self.freqCutoffLabel.place(x=20,y=220,width=80,height=20)

        #self.freqCutoffEntry = tk.Entry(top,justify='right')
        #self.freqCutoffEntry.place(x=20,y=240,width=100,height=30)
        #self.freqCutoffEntry.insert(0,"0.1")

        self.numThreadsLabel = tk.Label(top)
        self.numThreadsLabel.configure(text="Num. threads")
        self.numThreadsLabel.place(x=20,y=220,width=90,height=20)

        self.numThreadsEntry = tk.Entry(top,justify='right')
        self.numThreadsEntry.place(x=20,y=240,width=100,height=30)
        self.numThreadsEntry.insert(0,"8")
        
        self.runButton = tk.Button(top,command=mainAlgorithm)
        self.runButton.place(x=780,y=240,width=100,height=30)
        self.runButton.configure(text="Run")

        #self.plotButton = tk.Button(top)
        #self.plotButton.place(x=660,y=240,width=100,height=30)
        #self.plotButton.configure(text="Plot")

        self.exitButton = tk.Button(top,command=exitProgram)
        self.exitButton.place(x=660,y=240,width=100,height=30)
        self.exitButton.configure(text="Exit")
   
        self.logFileLabel = tk.Label(top)
        self.logFileLabel.place(x=20,y=280,height=20,width=100)
        self.logFileLabel.configure(text="Log window")

        self.logFrame = tk.Frame(top)
        self.logFrame.place(x=20, y=280, height=180, width=860)
        self.logFrame.configure(relief='groove')
        self.logFrame.configure(borderwidth="2")
        self.logFrame.configure(relief='groove')
        self.logFrame.configure(width=125)
        self.logArea = tk.Text(top,state='disabled')
        self.logArea.place(x=25,y=285,height=170, width=850)
        self.logArea.configure(background="white",borderwidth=5)
        self.logArea.configure(selectbackground="#c4c4c4")




if __name__ == '__main__':
    vp_start_gui()
