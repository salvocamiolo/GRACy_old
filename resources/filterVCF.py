#!/usr/bin/python
import sys

filename = sys.argv[1]

infile = open(filename)
outfile = open(filename+"_filtered.vcf","w")

line = infile.readline().rstrip()
outfile.write(line+"\n")

while line[0] == "#":
    line = infile.readline().rstrip()
    if not line:
        break
    if line[0]== "#":
        outfile.write(line+"\n")

while True:
    if not line:
        break

    allBeses= ((line.split('\t'))[4]).split(",")

    dp = float((((line.split("DP="))[1]).split(";"))[0])
    fields = line.split(":")
    alleleCounts = fields[-1]
    print alleleCounts
    alleles = alleleCounts.split(",")
    maxAllele = 0
    bestAllele = 0
    for a in range(0,len(alleles),+2):
        if (int(alleles[a]) + int(alleles[a+1]) ) > maxAllele:
            maxAllele = (int(alleles[a]) + int(alleles[a+1]) )
            bestAllele = a

    print "Best Allele is",maxAllele,bestAllele,allBeses[int(bestAllele/2)-1]
    if not bestAllele==0:
        split1 = line.split("\t")
        for a in range(4):
            outfile.write(split1[a]+"\t")
        outfile.write(allBeses[int(bestAllele/2)-1]+"\t")
        for a in range(5,len(split1)):
            outfile.write(split1[a]+"\t")
        outfile.write("\n")




   
    
    line = infile.readline().rstrip()

outfile.close()


