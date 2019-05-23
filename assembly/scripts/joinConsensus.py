from Bio import SeqIO
import sys

finalSequence = ""

for seq_record in SeqIO.parse("finalScaffold_1_15001_f.txt_con.fasta","fasta"):
    finalSequence+=str(seq_record.seq)

for seq_record in SeqIO.parse("finalScaffold_15001_200001_f.txt_con.fasta","fasta"):
    finalSequence+=str(seq_record.seq)

for seq_record in SeqIO.parse("finalScaffold_200001_2000000_f.txt_con.fasta","fasta"):
    finalSequence+=str(seq_record.seq)

outfile = open(sys.argv[1]+"_genome.fasta","w")
outfile.write(">finalScaffold\n"+finalSequence+"\n")

