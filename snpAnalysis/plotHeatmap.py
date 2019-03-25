#import seaborn as sns; sns.set()
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

position2PlotFile = sys.argv[1]
file2AnalyzeFile = sys.argv[2]

pos2plot = open(position2PlotFile)

outfile = open("snpsTable.txt","w")
outfile.write("Snp\t")

#Read the positions to plot
p2p = []
while True:
    line = pos2plot.readline().rstrip()
    if not line:
        break
    p2p.append(int(line))

p2p_sort = sorted(p2p)

heatmapFiles = open(file2AnalyzeFile)
snpPositions = {}
xlab = []
while True:
    line = heatmapFiles.readline().rstrip()
    if not line:
        break
    xlab.append(line)
    outfile.write(line+"\t")
    if not line in snpPositions:
        snpPositions[line] = {}
    infile = open(line)
    line2 = infile.readline().rstrip()
    while True:
        line2 = infile.readline().rstrip()
        if not line2:
            break
        fields = line2.split("\t")
        if not (fields[0],fields[1],fields[2]) in snpPositions[line]:
            snpPositions[line][(fields[0],fields[1],fields[2])] = float(fields[4])

    infile.close()
outfile.write("\n")

#Calculate the total number of snps present in the experiments
allSNPs = set()

for experiments in snpPositions:
    print "In the experiment",experiments,"There are",len(snpPositions[experiments]),"snps"
    for snp in snpPositions[experiments]:
        allSNPs.add(snp)

snp2Plot = []
for position in p2p_sort:
    for item in allSNPs:
        if int(item[0]) == position:
            snp2Plot.append(item)


#Create the table to plot
table2Plot = []
ylab=[]
for snp in snp2Plot:
    ylab.append(snp)
    outfile.write(str(snp)+"\t")
    line = []
    for experiments in snpPositions:
        if snp in snpPositions[experiments]:
            line.append(snpPositions[experiments][snp])
            outfile.write(str(snpPositions[experiments][snp])+"\t")
        else:
            line.append(np.nan)
            outfile.write("\t")
    table2Plot.append(line)
    outfile.write("\n")

print table2Plot
print len(table2Plot)
print xlab

#df = pd.DataFrame(table2Plot)
#df.columns = xlab
#print df
#fig,ax = 
#plt.subplots(figsize=(12,40))
#title = "Heatmap"
#plt.title(title,fontsize=18)
#ttl =ax.title
#ttl.set_position([0.5,1.05])
#ax.axis('off')

#sns.heatmap(table2Plot,xticklabels=xlab,yticklabels=ylab,annot_kws={"size":0.01})
#plt.show()
#plt.savefig('heatmap.png')
    






