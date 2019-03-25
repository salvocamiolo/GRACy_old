from Bio import SeqIO

kmerSeqs = {}

print "Load data in memory"
for seq_record in SeqIO.parse("counts.txt","fasta"):
    numOcc = int(str(seq_record.id))
    kmer = str(seq_record.seq)
    if not kmer in kmerSeqs:
        kmerSeqs[kmer] = numOcc

print "Calculate how many kmer has occurrence 1"
num1 = 0
for kmer in kmerSeqs:
    if kmerSeqs[kmer]==1:
        num1+=1

print "Result:", num1

