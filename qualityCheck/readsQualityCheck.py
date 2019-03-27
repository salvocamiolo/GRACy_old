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
from os import listdir
from os.path import isfile, join
import os
import time
import numpy as nm
import matplotlib.pyplot as plt


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


        def exitProgram():
           exit()

        def openInputFolder():
            inputFolder = tkFileDialog.askdirectory(initialdir = "~",title = "Select directory")
            self.inputFolderEntry.delete(0,tk.END)
            self.inputFolderEntry.insert(0,inputFolder)

        def openOutputFolder():
            outputolder = tkFileDialog.askdirectory(initialdir = "~",title = "Select directory")
            self.outputFolderEntry.delete(0,tk.END)
            self.outputFolderEntry.insert(0,outputolder)

        def openRecodingFile():
            recodingFile = tkFileDialog.askopenfilename(initialdir = "~",title = "Select file")
            self.recodingFileEntry.delete(0,tk.END)
            self.recodingFileEntry.insert(0,recodingFile)


        def plotCoveragePlot(covFile,expName):
            position = []
            coverage = []
            windowSize = 0
            avPos = 0
            covFile = open(covFile)
            coverageToPlot = []
            positionToPlot = []
            while True:
                line = covFile.readline().rstrip()
                if not line:
                    break
                fields = line.split("\t")
                position.append(int(fields[1])) 
                coverage.append(int(fields[2])) 
                windowSize +=1
                if windowSize == 200:
                    avPos += 200
                    positionToPlot.append(avPos)
                    avCoverage = nm.mean(coverage)
                    coverageToPlot.append(int(avCoverage))
                    position = []
                    coverage = []
                    windowSize = 0

            plt.figure(figsize=(20,10))
            plt.plot(positionToPlot,coverageToPlot,linewidth=0.3)
            plt.xticks(nm.arange(0,235000,10000))
            plt.xlabel("Position (bp)")
            plt.ylabel("Coverage")
            plt.text(1000,950,expName)
            plt.savefig(expName+"_covPlot.png")









#*******************************************************
#********************* Main algorithm start ************
#*******************************************************

        def runQualityCheck():
            inputFolder = self.inputFolderEntry.get()
            outputFolder = self.outputFolderEntry.get()

            #Collect files to quality check
            onlyfiles = [f for f in listdir(inputFolder) if isfile(join(inputFolder, f))]
            fileNumber = {}
            for item in onlyfiles:
                if not item[:-8] in fileNumber:
                    fileNumber[item[:-8]] = 1
                else:
                    fileNumber[item[:-8]] += 1

            paired2filter = []
            single2filter = []
            for files in fileNumber:
                if fileNumber[files]==2:
                    paired2filter.append(files)
                if fileNumber[files]==1:
                    single2filter.append(files)

            #Check recoding file
            codes = {}
            if self.recodeSampleChkValue.get() == True:
                rcodefile = open(self.recodingFileEntry.get())
                

                while True:
                    line = rcodefile.readline().rstrip()
                    if not line:
                        break
                    fields = line.split("\t")
                    oldName = (fields[0].split("/"))[-1]
                    if not oldName in codes:
                        codes[oldName]=fields[1]
            else:
                for item in paired2filter:
                    if not item+"_1.fastq" in codes:
                        codes[item+"_1.fastq"] = item+"_1.fastq"
                    if not item+"_2.fastq" in codes:
                        codes[item+"_2.fastq"] = item+"_2.fastq"


                    

            #Start filtering
            
            statisticsToReport = []
            statisticsToReport.append("Original reads")
            datasetStatistics = {}
            numTasks = 1
            if self.removeHumanChkValue.get() == True:
                numTasks += 4
                statisticsToReport.append("Number Reads after host reads removal")
                statisticsToReport.append("Percentage Reads after host reads removal")
            if self.removeAdaptersChkValue.get() == True:
                numTasks += 2
                statisticsToReport.append("Number Reads after adaptor trimming")
                statisticsToReport.append("Percentage Reads after adaptor trimming")
            if self.deduplicateChkValue.get() == True:
                numTasks += 2
                statisticsToReport.append("Number Reads after deduplication")
                statisticsToReport.append("Percentage Reads after deduplication")
            if self.al2RefChkValue.get() == True:
                numTasks += 5
                statisticsToReport.append("Merlin coverage for original reads")
                statisticsToReport.append("Number Original reads mapping to Merlin")
                statisticsToReport.append("Percentage Original reads mapping to Merlin")
                statisticsToReport.append("Number of bases covered by original reads")
                statisticsToReport.append("Percentage of bases covered by original reads")


                if self.deduplicateChkValue.get() == True:
                    numTasks += 5
                    statisticsToReport.append("Merlin coverage for deduplicated reads")
                    statisticsToReport.append("Number of deduplicated reads mapping to Merlin")
                    statisticsToReport.append("Percentage of deduplicated reads mapping to Merlin")
                    statisticsToReport.append("Number of bases covered by deduplicated reads")
                    statisticsToReport.append("Percentage of bases covered by deduplicated reads")




            numTasks = numTasks*len(paired2filter)
            progressBarIncrement = 100000/numTasks
            step = 0

            for dataset in paired2filter:
                if not dataset[0] ==".":
                    if not dataset in datasetStatistics:
                        datasetStatistics[dataset] = []
                    os.system("ln -s "+inputFolder+"/"+dataset+"_1.fastq tempReads_1.fastq")
                    os.system("ln -s "+inputFolder+"/"+dataset+"_2.fastq tempReads_2.fastq")

                    self.logArea.configure(state='normal')
                    self.logArea.insert(tk.END, "Calculating the number of reads for dataset "+dataset+"....\n")
                    self.logArea.see(tk.END)
                    self.logArea.configure(state='disabled')
                    self.logArea.update()
                    os.system("wc -l "+inputFolder+"/"+dataset+"_1.fastq  >numReads")
                    step+=1
                    self.progressbar['value']= int( (step*progressBarIncrement) )
                    numreads = open("numReads")
                    
                    numberOfReads = int(((numreads.readline().rstrip()).split(" "))[0])/2
                    originalReadsNumber = numberOfReads
                    datasetStatistics[dataset].append(numberOfReads)
                    numreads.close()
                    os.system("rm -f numReads")
                    self.logArea.configure(state='normal')
                    self.logArea.insert(tk.END, "Done!...... Number of reads: "+str(numberOfReads)+"\n\n")
                    self.logArea.see(tk.END)
                    self.logArea.configure(state='disabled')
                    self.logArea.update()


                    suffixCode = ""
                    if self.removeHumanChkValue.get() == True:
                        suffixCode+="_nh"
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Removing human reads from dataset "+dataset+"\n")
                        self.logArea.configure(state='disabled')
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Mapping reads "+dataset+" to the host reference genome....")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        os.system("bowtie2 --local -x /home2/db/bowtie2/hg38 -1 tempReads_1.fastq -2 tempReads_2.fastq -p "+self.numThreadsEntry.get()+" -S hostAlignment.sam")
                        step+=1
                        self.progressbar['value']= int( (step*progressBarIncrement) )
                        self.progressbar.update()
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Done!\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
    
                        step+=1
                        self.progressbar['value']= int( (step*progressBarIncrement) )
                        self.progressbar.update()
                        time.sleep(2)
                        self.progressbar.update()
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Converting alignment format for reads "+dataset+"....\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        os.system("samtools view -bS -h hostAlignment.sam >hostAlignment.bam") 
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Done!\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        step+=1
                        self.progressbar['value']= int( (step*progressBarIncrement) )
                        self.progressbar.update()
                        time.sleep(2)
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Extracting unmapped reads for dataset "+dataset+"....")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        os.system("bam2fastq --no-aligned --force --strict -o unmapped#.fq hostAlignment.bam")
                        step+=1
                        self.progressbar['value']= int( (step*progressBarIncrement) )
                        self.progressbar.update()
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Done!\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        os.system("mv unmapped_1.fq tempReads_1.fastq")
                        os.system("mv unmapped_2.fq tempReads_2.fastq")
                        os.system("rm -f hostAlignment.bam ")

                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Calculating the number of host free reads for dataset "+dataset+"....")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        os.system("wc -l tempReads_1.fastq >numReads")
                        step+=1
                        self.progressbar['value']= int( (step*progressBarIncrement) )
                        self.progressbar.update()
                        numreads = open("numReads")
                        numberOfReads = int(((numreads.readline().rstrip()).split(" "))[0])/2
                        datasetStatistics[dataset].append(str(numberOfReads))
                        datasetStatistics[dataset].append(str(numberOfReads*100/originalReadsNumber)[:4])
                        numreads.close()
                        os.system("rm -f numReads")
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Done!...... Number of reads: "+str(numberOfReads)+"\n\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
    


                        #step += 1
                        #self.progressbar['value']=int( (step/numStep)*100)
                        self.progressbar.update()
                        time.sleep(2)     

                    if self.removeAdaptersChkValue.get() == True:
                        suffixCode += "_tr"
        
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Trimming reads for dataset "+dataset+"....")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        print "Performing trim_galore"
                        if len(self.adapter1Entry.get())>5 and len(self.adapter2Entry.get())>5:
                            adpt1 = self.adapter1Entry.get()
                            adpt2 = self.adapter2Entry.get()
                            os.system("trim_galore -paired -a "+adpt1+" -a2 "+adpt2+" tempReads_1.fastq tempReads_2.fastq ")

                        else:
                            os.system("trim_galore -paired tempReads_1.fastq tempReads_2.fastq ")

                        step+=1
                        self.progressbar['value']= int( (step*progressBarIncrement) )
                        self.progressbar.update()
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Done!\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()

                        os.system("mv tempReads_1_val_1.fq tempReads_1.fastq")
                        os.system("mv tempReads_2_val_2.fq tempReads_2.fastq")
                        os.system("rm -f tempReads_1.fastq_trimming_report.txt tempReads_2.fastq_trimming_report.txt")

                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Calculating the number reads after trimming for dataset "+dataset+"....")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        os.system("wc -l tempReads_1.fastq >numReads")
                        step+=1
                        self.progressbar['value']= int( (step*progressBarIncrement) )
                        self.progressbar.update()
                        numreads = open("numReads")
                        numberOfReads = int(((numreads.readline().rstrip()).split(" "))[0])/2
                        datasetStatistics[dataset].append(str(numberOfReads))
                        datasetStatistics[dataset].append(str(numberOfReads*100/originalReadsNumber)[:4])
                        numreads.close()
                        os.system("rm -f numReads")
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Done!...... Number of reads: "+str(numberOfReads)+"\n\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        



                    if self.deduplicateChkValue.get() == True:

                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Performing deduplication on dataset "+dataset+"....")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        os.system("fastuniq -i dedupTemplate.txt -t q -o "+codes[dataset+"_1.fastq"].replace("_1.fastq","")+suffixCode+"_dd_1.fastq -p "+codes[dataset+"_2.fastq"].replace("_2.fastq","")+suffixCode+"_dd_2.fastq")
                        dedupFileName1 = codes[dataset+"_1.fastq"].replace("_1.fastq","")+suffixCode+"_dd_1.fastq"
                        dedupFileName2 = codes[dataset+"_2.fastq"].replace("_2.fastq","")+suffixCode+"_dd_2.fastq"
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Done!\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        step+=1
                        self.progressbar['value']= int( (step*progressBarIncrement) )
                        self.progressbar.update()
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Calculating the number reads after deduplicating for dataset "+dataset+"....")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        os.system("wc -l "+dedupFileName1+" >numReads")
                        step+=1
                        self.progressbar['value']= int( (step*progressBarIncrement) )
                        self.progressbar.update()
                        numreads = open("numReads")
                        numberOfReads = int(((numreads.readline().rstrip()).split(" "))[0])/2
                        datasetStatistics[dataset].append(str(numberOfReads))
                        datasetStatistics[dataset].append(str(numberOfReads*100/originalReadsNumber)[:4])
                        numreads.close()
                        os.system("rm -f numReads")
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Done!...... Number of reads: "+str(numberOfReads)+"\n\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        

                    if self.al2RefChkValue.get() == True:
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Performing alignment for original dataset "+dataset+"....")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        os.system("bowtie2 -1 tempReads_1.fastq -2 tempReads_2.fastq -x ../referenceSequences/hcmv -p "+self.numThreadsEntry.get()+" -S alignment.sam")
                        step+=1
                        self.progressbar['value']= int( (step*progressBarIncrement) )
                        self.progressbar.update()
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Done!\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        
                        
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Convert sam to bam....")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        os.system("samtools view -bS -h  alignment.sam >alignment.bam")
                        step+=1
                        self.progressbar['value']= int( (step*progressBarIncrement) )
                        self.progressbar.update()
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Done!\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()


                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Sorting bam....")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        os.system("samtools sort -o alignment_sorted.bam alignment.bam")
                        step+=1
                        self.progressbar['value']= int( (step*progressBarIncrement) )
                        self.progressbar.update()
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Done!\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()


                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Calculating coverage for original dataset "+dataset+"....")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        os.system("samtools depth -d 10000000 alignment_sorted.bam  |  awk '{sum+=$3} END { print sum/NR}' >avCoverage.txt")
                        os.system("samtools depth -a -d 10000000 alignment_sorted.bam >coverage.txt")
                        plotCoveragePlot("coverage.txt",codes[dataset+"_1.fastq"].replace("_1.fastq","")+"_nh_tr")
                        os.system("samtools depth -d 10000000 alignment_sorted.bam  | wc -l >breadth.txt")
                        step+=1
                        self.progressbar['value']= int( (step*progressBarIncrement) )
                        self.progressbar.update()
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Done!\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()


                        avCovFile = open("avCoverage.txt")
                        avCov = float(avCovFile.readline().rstrip())
                        avCovFile.close()
                        print "Average coverage",avCov
                        
                        datasetStatistics[dataset].append(avCov)

                        os.system("samtools view -F 4 alignment.sam | wc -l > readsMapping")
                        step+=1
                        self.progressbar['value']= int( (step*progressBarIncrement) )
                        self.progressbar.update()
                        readsMappingFile = open("readsMapping")
                        readsMapping = float(readsMappingFile.readline().rstrip())
                        readsMappingFile.close()
                        
                        datasetStatistics[dataset].append(str(readsMapping))
                        datasetStatistics[dataset].append(str(readsMapping*100/originalReadsNumber)[:4])

                        bfile = open("breadth.txt")
                        breadthValue = float(bfile.readline().rstrip())
                        bfile.close()

                        datasetStatistics[dataset].append(str(breadthValue))
                        datasetStatistics[dataset].append(str( breadthValue / 235646.0 )[:4])



                        if self.deduplicateChkValue.get() == True:
                            self.logArea.configure(state='normal')
                            self.logArea.insert(tk.END, "Performing alignment for deduplicated dataset "+dataset+"....")
                            self.logArea.see(tk.END)
                            self.logArea.configure(state='disabled')
                            self.logArea.update()
                            os.system("bowtie2 -1 "+dedupFileName1+" -2 "+ dedupFileName2+ "  -x ../referenceSequences/hcmv -p "+self.numThreadsEntry.get()+" -S alignment.sam")
                            step+=1
                            self.progressbar['value']= int( (step*progressBarIncrement) )
                            self.progressbar.update()
                            self.logArea.configure(state='normal')
                            self.logArea.insert(tk.END, "Done!\n")
                            self.logArea.see(tk.END)
                            self.logArea.configure(state='disabled')
                            self.logArea.update()
                            
                            self.logArea.configure(state='normal')
                            self.logArea.insert(tk.END, "Converting sam to bam....")
                            self.logArea.see(tk.END)
                            self.logArea.configure(state='disabled')
                            self.logArea.update()
                            os.system("samtools view -bS -h  alignment.sam >alignment.bam")
                            step+=1
                            self.progressbar['value']= int( (step*progressBarIncrement) )
                            self.progressbar.update()
                            self.logArea.configure(state='normal')
                            self.logArea.insert(tk.END, "Done!\n")
                            self.logArea.see(tk.END)
                            self.logArea.configure(state='disabled')
                            self.logArea.update()
                            
                            self.logArea.configure(state='normal')
                            self.logArea.insert(tk.END, "Sorting bam file....")
                            self.logArea.see(tk.END)
                            self.logArea.configure(state='disabled')
                            self.logArea.update()
                            os.system("samtools sort -o alignment_sorted.bam alignment.bam")
                            step+=1
                            self.progressbar['value']= int( (step*progressBarIncrement) )
                            self.progressbar.update()
                            self.logArea.configure(state='normal')
                            self.logArea.insert(tk.END, "Done!\n")
                            self.logArea.see(tk.END)
                            self.logArea.configure(state='disabled')
                            self.logArea.update()
                            
                            self.logArea.configure(state='normal')
                            self.logArea.insert(tk.END, "Calculating coverage for deduplicated dataset "+dataset+"....")
                            self.logArea.see(tk.END)
                            self.logArea.configure(state='disabled')
                            self.logArea.update()
                            os.system("samtools depth -d 10000000 alignment_sorted.bam  |  awk '{sum+=$3} END { print sum/NR}' >avCoverage.txt")
                            os.system("samtools depth -a -d 10000000 alignment_sorted.bam   >coverage.txt")
                            plotCoveragePlot("coverage.txt",codes[dataset+"_1.fastq"].replace("_1.fastq","")+"_nh_tr_dd")
                            os.system("samtools depth -d 10000000 alignment_sorted.bam  | wc -l >breadth.txt")
                            step+=1
                            self.progressbar['value']= int( (step*progressBarIncrement) )
                            self.progressbar.update()
                            avCovFile = open("avCoverage.txt")
                            avCov = float(avCovFile.readline().rstrip())
                            avCovFile.close()
                            self.logArea.configure(state='normal')
                            self.logArea.insert(tk.END, "Done!\n")
                            self.logArea.see(tk.END)
                            self.logArea.configure(state='disabled')
                            self.logArea.update()

                            print "Average coverage for deduplicated ",avCov
                            
                            datasetStatistics[dataset].append(avCov)

                            os.system("samtools view -F 4 alignment.sam | wc -l > readsMapping")
                            step+=1
                            self.progressbar['value']= int( (step*progressBarIncrement) )
                            self.progressbar.update()
                            readsMappingFile = open("readsMapping")
                            readsMapping = float(readsMappingFile.readline().rstrip())
                            readsMappingFile.close()
                            
                            datasetStatistics[dataset].append(str(readsMapping))
                            datasetStatistics[dataset].append(str(readsMapping*100/originalReadsNumber)[:4])

                            bfile = open("breadth.txt")
                            breadthValue = float(bfile.readline().rstrip())
                            bfile.close()

                            datasetStatistics[dataset].append(str(breadthValue))
                            datasetStatistics[dataset].append(str( breadthValue / 235646.0 )[:4])





                    print "finished"
                    os.system("mv tempReads_1.fastq "+outputFolder+"/"+codes[dataset+"_1.fastq"].replace("_1.fastq","")+suffixCode+"_1.fastq")
                    os.system("mv tempReads_2.fastq "+outputFolder+"/"+codes[dataset+"_2.fastq"].replace("_2.fastq","")+suffixCode+"_2.fastq")
                    if self.deduplicateChkValue.get() == True:
                        os.system("mv "+dedupFileName1+" "+outputFolder+"/")
                        os.system("mv "+dedupFileName2+" "+outputFolder+"/")
                    os.system("mv "+codes[dataset+"_1.fastq"].replace("_1.fastq","")+"* "+outputFolder+"/")
                    
                    for a in range(len(statisticsToReport)):
                        print statisticsToReport[a],datasetStatistics[dataset][a]

                    os.system("rm -f *.sam *.bam readsMapping avCoverage.txt breadth.txt")



                    for item in datasetStatistics:
                        print item
                        for a in datasetStatistics[item]:
                            print str(a)+"\t"
                    print "\n" 


            outfile = open("summaryTable.txt","w")
            outfile.write("Statistics\t")


            for item in datasetStatistics:
                outfile.write(codes[item+"_1.fastq"].replace("_1.fastq","")+"\t")
            outfile.write("\n")
            for a in range(len(statisticsToReport)):
                outfile.write(statisticsToReport[a]+"\t")
                for item in datasetStatistics:
                    outfile.write(str(datasetStatistics[item][a])+"\t")
                outfile.write("\n")
                

               


        


        '''This class configures and populates the selflevel window.
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

        w=880
        h=660
        ws = root.winfo_screenwidth() # width of the screen
        hs = root.winfo_screenheight() # height of the screen
        # calculate x and y coordinates for the Tk root window
        x = (ws/2) - (w/2)
        y = (hs/2) - (h/2)

        top.geometry('%dx%d+%d+%d' % (w, h, x, y)) 
        top.title("Quality checker")
        top.configure(highlightcolor="black")


        

        self.inputFolderLabel = tk.Label(top)
        self.inputFolderLabel.place(x=20,y=20,height=20, width=80)
        self.inputFolderLabel.configure(text="Input folder")

        self.inputFolderEntry = tk.Entry(top)
        self.inputFolderEntry.place(x=20,y=40,height=30, width = 300)
        #self.inputFolderEntry.insert(0,"Please select folder....")
        self.inputFolderEntry.insert(0,"/home3/scc20x/infolder1")

        self.inputFolderButton = tk.Button(top,command=openInputFolder)
        self.inputFolderButton.place(x=340,y=40,height=30,width=100)
        self.inputFolderButton.configure(text="Open folder")

        self.outputFolderLabel = tk.Label(top)
        self.outputFolderLabel.place(x=20,y=80,height=20, width=90)
        self.outputFolderLabel.configure(text="Output folder")

        self.outputFolderEntry = tk.Entry(top)
        self.outputFolderEntry.place(x=20,y=100,height=30, width = 300)
        #self.outputFolderEntry.insert(0,"Please select folder....")
        self.outputFolderEntry.insert(0,"/home3/scc20x/outfolder1")

        self.outputFolderButton = tk.Button(top,command=openOutputFolder)
        self.outputFolderButton.place(x=340,y=100,height=30,width=100)
        self.outputFolderButton.configure(text="Open folder")

        self.recodingFileLabel = tk.Label(top,justify="left")
        self.recodingFileLabel.place(x=20,y=140,height=20,width=90)
        self.recodingFileLabel.configure(text="Recoding file")

        self.recodingFileEntry = tk.Entry(top)
        self.recodingFileEntry.place(x=20,y=160,height=30, width = 300)
        #self.recodingFileEntry.insert(0,"Please select file....")
        self.recodingFileEntry.insert(0,"./recodingFile.txt")



        self.recodingFileButton = tk.Button(top,command=openRecodingFile)
        self.recodingFileButton.place(x=340,y=160,height=30,width=100)
        self.recodingFileButton.configure(text="Open file")

        self.adapter1Label = tk.Label(top)
        self.adapter1Label.place(x=20,y=200,height=20, width=70)
        self.adapter1Label.configure(text="Adapter 1")

        self.adapter1Entry = tk.Entry(top)
        self.adapter1Entry.place(x=20,y=220,height=30, width = 420)

        self.adapter2Label = tk.Label(top)
        self.adapter2Label.place(x=20,y=260,height=20, width=70)
        self.adapter2Label.configure(text="Adapter 2")

        self.adapter2Entry = tk.Entry(top)
        self.adapter2Entry.place(x=20,y=280,height=30, width = 420)


        self.numThreadsLabel = tk.Label(top,justify="left")
        self.numThreadsLabel.place(x=20,y=320,height=20,width=130)
        self.numThreadsLabel.configure(text="Number of threads")

        self.numThreadsEntry = tk.Entry(top,justify="right")
        self.numThreadsEntry.place(x=20,y=340,height=30, width = 150)
        self.numThreadsEntry.insert(0,"8")

        self.tasksLabel = tk.Label(top)
        self.tasksLabel.place(x=500,y=20,height=20,width=120)
        self.tasksLabel.configure(text="Tasks to perform")

        self.tasksFrame = tk.Frame(top)
        self.tasksFrame.place(x=500, y=40, height=220, width=340)
        self.tasksFrame.configure(borderwidth="2")
        self.tasksFrame.configure(relief='groove')

        self.removeHumanChkValue = tk.BooleanVar() 
        self.removeHumanChkValue.set(False)
        self.removeHumanCheckButton = tk.Checkbutton(top,variable=self.removeHumanChkValue)
        self.removeHumanCheckButton.place(x=520,y=60,height=20,width=170)
        self.removeHumanCheckButton.configure(text="Remove human reads")

        self.removeAdaptersChkValue = tk.BooleanVar() 
        self.removeAdaptersChkValue.set(False)
        self.removeAdaptersCheckButton = tk.Checkbutton(top,variable=self.removeAdaptersChkValue)
        self.removeAdaptersCheckButton.place(x=520,y=100,height=20,width=115)
        self.removeAdaptersCheckButton.configure(text="Trim adapters")

        self.recodeSampleChkValue = tk.BooleanVar() 
        self.recodeSampleChkValue.set(False)
        self.recodeSamplesCheckButton = tk.Checkbutton(top,variable=self.recodeSampleChkValue)
        self.recodeSamplesCheckButton.place(x=520,y=140,height=20,width=170)
        self.recodeSamplesCheckButton.configure(text="Recode sample names")

        self.deduplicateChkValue = tk.BooleanVar() 
        self.deduplicateChkValue.set(False)
        self.deduplicateCheckButton = tk.Checkbutton(top,variable=self.deduplicateChkValue)
        self.deduplicateCheckButton.place(x=520,y=180,height=20,width=220)
        self.deduplicateCheckButton.configure(text="Deduplicate clonal fragments")

        self.al2RefChkValue = tk.BooleanVar() 
        self.al2RefChkValue.set(False)
        self.al2RefCheckButton = tk.Checkbutton(top,variable=self.al2RefChkValue)
        self.al2RefCheckButton.place(x=520,y=220,height=20,width=190)
        self.al2RefCheckButton.configure(text="Align to Merlin reference")


        self.logTextLabel = tk.Label(top)
        self.logTextLabel.place(x=20,y=400,height=20,width=90)
        self.logTextLabel.configure(text="Log window")
        self.logFrame = tk.Frame(top)
        self.logFrame.place(x=20, y=420, height=180, width=820)
        self.logFrame.configure(relief='groove')
        self.logFrame.configure(borderwidth="2")
        self.logFrame.configure(relief='groove')
        self.logFrame.configure(width=125)
        self.logArea = tk.Text(top,state='disabled')
        self.logArea.place(x=25,y=425,height=170, width=810)
        self.logArea.configure(background="white",borderwidth=5)
        self.logArea.configure(selectbackground="#c4c4c4")

        self.progressbar=ttk.Progressbar(top,orient="horizontal",length=820,mode="determinate")
        self.progressbar.place(x=20,y=620)
        self.progressbar['maximum'] = 100000

        self.runButton = tk.Button(top,command=runQualityCheck)
        self.runButton.place(x=740,y=360,height=30,width=100)
        self.runButton.configure(text="Run")

        self.exitButton = tk.Button(top,command=exitProgram)
        self.exitButton.place(x=620,y=360,height=30,width=100)
        self.exitButton.configure(text="Exit")









if __name__ == '__main__':
    vp_start_gui()
