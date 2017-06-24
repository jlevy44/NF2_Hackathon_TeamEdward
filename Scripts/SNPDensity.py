from pybedtools import BedTool
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

with open('/mnt/disks/data-bam-tumor/Homo_sapiens_assembly19.fasta.fai','r') as f:
    chromosomes = {line.split('\t')[0]:int(line.split('\t')[1]) for line in f.readlines()[0:24]}

print chromosomes

chrom_Intervals = {k:np.append(np.arange(0,chromosomes[k],10000),chromosomes[k]) for k in chromosomes.keys()}

for k in chromosomes.keys():
    print chrom_Intervals[k][0:10]


with open('humanBedInterval.bed','w') as f:
    for k in chromosomes.keys():
        f.write('\n'.join('%s\t%d\t%d'%(k,chrom_Intervals[k][i],chrom_Intervals[k][i+1]-1) for i in range(len(chrom_Intervals[k])-1)) +'\n')

a = BedTool('humanBedInterval.bed').sort()
b = BedTool('/mnt/disks/data-vcf/GSN79Tumor_normal.vcf')
print a.head()
print b.head()

a.coverage(b).saveas('VCFCoverage.bed')#,hist=True

#for k in chromosomes.keys():
positionHistogram = defaultdict(list)

with open('/mnt/disks/data-vcf/GSN79Tumor_normal.vcf','r') as f:
    for line in f.readlines():
        if line and line.startswith('#') ==0:
            if line.split('\t')[0] not in positionHistogram.keys():
                positionHistogram[line.split('\t')[0]] = []
            positionHistogram[line.split('\t')[0]].append(int(line.split('\t')[1]))

for k in chromosomes.keys():
    plt.hist(positionHistogram[k],bins=chrom_Intervals[k])
    plt.savefig(k+'SNPDense.png')