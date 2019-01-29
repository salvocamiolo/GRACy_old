from Bio import SeqIO
import biomodule
import os
import sys

read1 = sys.argv[1]
read2 = sys.argv[2]

toAssemble = []

exclude = 1
if exclude == 0:
    #Analyze scaffold builder output file and retrieve corresponding scaffolds
    infile = open("sb_output.txt")
    line = infile.readline()
    while not line[:7] == "#contig":
        line = infile.readline()

    while True:
        line = infile.readline().rstrip()
        if not line:
            break
        split1 = line.split("\t")
        split2 = []
        for item in split1:
            if not item == '':
                split2.append(item)
        title = split2[0]
        start = int(split2[3])
        end = int(split2[4])
        if int(split2[5]) > 500:
            toAssemble.append(title+".fasta")
            os.system("./getSequenceFromFasta.py ../2_spadesAssembly/scaffolds.fasta "+title+" "+title+".fasta")
            if start > end:
                os.system("revcomp "+title+".fasta > temp.fasta; mv temp.fasta "+title+".fasta")

    #print toAssemble

    #Join scaffolds

    #If there is more than one contig they will all be joined together and the final scaffold called finalScaffold_1.fasta
    #This sequence will be joined with the sequence_a.fasta in the reverse way and produce a file called finalScaffold_2.fasta
    #This sequence will be joined to the final sequence_a in forward and produce a file named finalScaffold.fasta
    if len(toAssemble)>1:
        for a in range(len(toAssemble)-1):
            while True:
                print "Start joining "+toAssemble[a]+" f "+toAssemble[a+1]+" f"
                #sys.stdin.read(1)
                print "Performing joinScaffold_careful algorithm...."
                os.system("python joinScaffolds_careful.py join "+read1+" "+read2+" "+toAssemble[a]+" f "+toAssemble[a+1]+" f")
                if os.path.isfile("joined_"+toAssemble[a]+"_"+toAssemble[a+1]) == True:
                    os.system("cp joined_"+toAssemble[a]+"_"+toAssemble[a+1]+" "+toAssemble[a+1])
                    break
                else:
                    print "WARNING! The joinScaffold algorithm did not work. Trying joinScaffold"
                    os.system("python joinScaffolds.py join "+read1+" "+read2+" "+toAssemble[a]+" f "+toAssemble[a+1]+" f")
                    if os.path.isfile("joined_"+toAssemble[a]+"_"+toAssemble[a+1]) == True:
                        os.system("cp joined_"+toAssemble[a]+"_"+toAssemble[a+1]+" "+toAssemble[a+1])
                        break
                    else:
                        print "Also joinScaffold_careful failed! Now trying joinScaffold_trivial: 10 steps....."
                        inputFile = "joined_W_WARNING_"+toAssemble[a]+"_"+toAssemble[a+1]
                        for seq_record in SeqIO.parse(inputFile, "fasta"):
                            inputSeq = str(seq_record.seq)
                            initialLen = len(inputSeq)
                        print "Initial length",initialLen
                        os.system("python joinScaffolds_trivial.py "+read1+" "+read2+" "+inputFile+" f " + toAssemble[a+1]+ " f 10 ")
                        for seq_record in SeqIO.parse("joinScaffold_trivialSeq.fasta", "fasta"):
                            outputSeq = str(seq_record.seq)
                            outputLen = len(outputSeq)
                        print "Final length",outputLen

                        if outputLen == initialLen:
                            print "RANGE NOT CLOSED ",toAssemble[a],toAssemble[a+1]
                            os.system("touch criticalError")
                            exit()
                        else:
                            os.system("mv joinScaffold_trivialSeq.fasta "+toAssemble[a])





                
            
            
            
            print "End joining "+toAssemble[a]+" f "+toAssemble[a+1]+" f"
            #sys.stdin.read(1)

        os.system("cp joined_"+toAssemble[a]+"_"+toAssemble[a+1]+" finalScaffold_1.fasta")
        os.system("lastz ~/hcmvReference/hcmv_genome.fasta finalScaffold_1.fasta --format=general:name1,strand1,start1,end1,name2,strand2,start2,end2 >lastzResult.txt")

        
        #sys.stdin.read(1)
    else:
        os.system("cp "+toAssemble[0]+" finalScaffold_1.fasta")

    fs = open("temp.fasta","w")
    for seq_record in SeqIO.parse("finalScaffold_1.fasta","fasta"):
        fs.write(">finalScaffold\n"+(str(seq_record.seq))[200:-200]+"\n")
    fs.close()
    os.system("mv temp.fasta finalScaffold_1.fasta")

else:
    #finalScaffoldFile = open("finalScaffold_1.fasta","w")
    #sbFile = open("sb_Scaffold.fasta")
    #sbSequence = sbFile.readline()
    #sbSequence = sbFile.readline().rstrip()
    #finalScaffoldFile.write(">finalScaffold\n"+sbSequence+"\n")
    #finalScaffoldFile.close()






    #Attach sequence a at 5prime end
    for seq_record in SeqIO.parse("finalScaffold_1.fasta","fasta"):
        startSeq = (str(seq_record.seq))[:10000]
        sequence = str(seq_record.seq)

    startSeqFile = open("startSeq.fasta","w")
    startSeqFile.write(">startSeq\n"+startSeq+"\n")
    startSeqFile.close()
    os.system("makeblastdb -dbtype nucl -in startSeq.fasta")
    os.system("blastn -query sequence_a.fasta -db startSeq.fasta -task blastn -outfmt 6 -out outputBlast.txt")
    outblast = open("outputBlast.txt")
    line = outblast.readline().rstrip()
    blastRange = 0
    if line:
        fields=line.split("\t")
        blastRange = (int(fields[7]) - int(fields[6]))

    if blastRange <20 or int(fields[8]>2000):
        print "Searching for 5 prime end...."
        if os.path.isfile("file5prime.fasta"):
            os.system("cp file5prime.fasta sequence_a.fasta")
        numTry = 0
        while True:
            numTry+=1
            if numTry>1:
                print "Could not get the 5 prime end"
                os.system("cp finalScaffold_1.fasta finalScaffold_2.fasta")
                break
            print "Performing joinScaffold_careful...."
            os.system("python joinScaffolds_careful.py join "+read1+" "+read2+" startSeq.fasta r sequence_a.fasta r")
            if os.path.isfile("joined_startSeq.fasta_sequence_a.fasta") == True:
                for seq_record in SeqIO.parse("joined_startSeq.fasta_sequence_a.fasta","fasta"):
                    newBit = biomodule.reverseComplement(str(seq_record.seq))
                sequence = newBit+sequence[10000:]
                finalScaffoldFile = open("finalScaffold_2.fasta","w")
                finalScaffoldFile.write(">finalScaffold\n"+sequence+"\n")
                finalScaffoldFile.close()
                break
            else:
                print "JoinScaffold_careful did not work! Trying joinScaffold...."
                os.system("cp toElong.fasta bestJoinScaffoldCarefulGuess.fasta")
                os.system("python joinScaffolds.py join "+read1+" "+read2+" startSeq.fasta r sequence_a.fasta r")
                if os.path.isfile("joined_startSeq.fasta_sequence_a.fasta") == True:
                    for seq_record in SeqIO.parse("joined_startSeq.fasta_sequence_a.fasta","fasta"):
                        newBit = biomodule.reverseComplement(str(seq_record.seq))
                    sequence = newBit+sequence[10000:]
                    finalScaffoldFile = open("finalScaffold_2.fasta","w")
                    finalScaffoldFile.write(">finalScaffold\n"+sequence+"\n")
                    finalScaffoldFile.close()
                    break
                else:
                    print "Also joinScaffold_careful failed! Now trying joinScaffold_trivial: 20 steps....."
                    inputFile = "bestJoinScaffoldCarefulGuess.fasta"
                    
                    for seq_record in SeqIO.parse(inputFile, "fasta"):
                        inputSeq = str(seq_record.seq)
                        initialLen = len(inputSeq)
                    print "Initial length",initialLen
                    os.system("python joinScaffolds_trivial.py "+read1+" "+read2+" "+inputFile+" f sequence_a.fasta r 20")
                    for seq_record in SeqIO.parse("joinScaffold_trivialSeq.fasta", "fasta"):
                        outputSeq = str(seq_record.seq)
                        outputLen = len(outputSeq)
                    print "Final length",outputLen
                    os.system("revcomp joinScaffold_trivialSeq.fasta >temp; mv temp joinScaffold_trivialSeq.fasta")
                    os.system("mv joinScaffold_trivialSeq.fasta 5primeEnd_failedReconstruction")

                    #if outputLen == initialLen:
                    #    print "RANGE NOT CLOSED finalScaffold_1 sequence_a"
                    #    os.system("touch criticalError")
                    #    exit()
                    #else:
                    #    os.system("revcomp joinScaffold_trivialSeq.fasta >temp; mv temp joinScaffold_trivialSeq.fasta")
                    #    os.system("mv joinScaffold_trivialSeq.fasta startSeq.fasta.fasta")





        #os.system("revcomp joined_finalScaffold_1.fasta_sequence_a.fasta > finalScaffold_2.fasta")

        print "End Attach 5 prime sequence A"
        #sys.stdin.read(1)
        #os.system("python joinScaffolds.py join "+read1+" "+read2+" finalScaffold_2.fasta f sequence_a.fasta f")
        #os.system("cp joined_finalScaffold_2.fasta_sequence_a.fasta finalScaffold.fasta")
    else:
        finalScaffoldFile = open("finalScaffold_2.fasta","w")
        finalScaffoldFile.write(">finalScaffold\n"+sequence[int(fields[8]):])
        finalScaffoldFile.close()

    #Attach sequence a at 3 prime end
    for seq_record in SeqIO.parse("finalScaffold_2.fasta","fasta"):
        startSeq = (str(seq_record.seq))[-10000:]
        sequence = str(seq_record.seq)

    startSeqFile = open("startSeq.fasta","w")
    startSeqFile.write(">startSeq\n"+startSeq+"\n")
    startSeqFile.close()

    os.system("makeblastdb -dbtype nucl -in startSeq.fasta")
    os.system("blastn -query sequence_a_prime.fasta -db startSeq.fasta -task blastn -outfmt 6 -out outputBlast.txt")
    outblast = open("outputBlast.txt")
    line = outblast.readline().rstrip()
    blastRange = 0
    if line:
        fields=line.split("\t")
        blastRange = (int(fields[7]) - int(fields[6]))
    #print blastRange
    #sys.stdin.read(1)


    if blastRange <20:
        print "Searching for 3 prime end...."
        if os.path.isfile("file3prime.fasta"):
            os.system("cp file3prime.fasta sequence_a_prime.fasta")
        
        numTry = 0
        while True:
            numTry += 1
            if numTry >1:
                print "Three prime end could not be reached...."
                os.system("cp finalScaffold_2.fasta finalScaffold.fasta")
                break
            print "Performing joinScaffold_careful...."
            os.system("python joinScaffolds_careful.py join "+read1+" "+read2+" startSeq.fasta f sequence_a_prime.fasta f")
            if os.path.isfile("joined_startSeq.fasta_sequence_a_prime.fasta") == True:
                for seq_record in SeqIO.parse("joined_startSeq.fasta_sequence_a_prime.fasta","fasta"):
                    newBit = str(seq_record.seq)
                sequence = sequence[:-10000]+newBit
                finalScaffoldFile = open("finalScaffold.fasta","w")
                finalScaffoldFile.write(">finalScaffold\n"+sequence+"\n")
                finalScaffoldFile.close()

                break
            else:
                print "JoinScaffold_careful did not work! Trying joinScaffold...."
                os.system("python joinScaffolds.py join "+read1+" "+read2+" startSeq.fasta f sequence_a_prime.fasta f")
                if os.path.isfile("joined_startSeq.fasta_sequence_a_prime.fasta") == True:
                    for seq_record in SeqIO.parse("joined_startSeq.fasta_sequence_a_prime.fasta","fasta"):
                        newBit = str(seq_record.seq)
                    sequence = sequence[:-10000]+newBit
                    finalScaffoldFile = open("finalScaffold.fasta","w")
                    finalScaffoldFile.write(">finalScaffold\n"+sequence+"\n")
                    finalScaffoldFile.close()
                    break
                else:
                    print "Also joinScaffold_careful failed! Now trying joinScaffold_trivial: 10 steps....."
                    inputFile = "joined_W_WARNING_startSeq.fasta_sequence_a_prime.fasta"
                    
                    for seq_record in SeqIO.parse(inputFile, "fasta"):
                        inputSeq = str(seq_record.seq)
                        initialLen = len(inputSeq)
                    print "Initial length",initialLen
                    os.system("python joinScaffolds_trivial.py "+read1+" "+read2+" "+inputFile+" f sequence_a_prime.fasta f 10")
                    for seq_record in SeqIO.parse("joinScaffold_trivialSeq.fasta", "fasta"):
                        outputSeq = str(seq_record.seq)
                        outputLen = len(outputSeq)
                    print "Final length",outputLen
                    os.system("mv joinScaffold_trivialSeq.fasta 3primeEnd_failedReconstruction")

                    #if outputLen == initialLen:
                    #    print "RANGE NOT CLOSED finalScaffold_2 sequence_a_prime"
                    #    os.system("touch criticalError")
                    #    exit()
                    #else:
                    #    os.system("mv joinScaffold_trivialSeq.fasta startSeq.fasta")

    else:
        finalScaffoldFile = open("finalScaffold.fasta","w")
        finalScaffoldFile.write(">finalScaffold\n"+sequence[:len(sequence)-10000+int(fields[9])])
        #print "Arrivato"
        finalScaffoldFile.close()
        #sys.stdin.read(1)
        





 





