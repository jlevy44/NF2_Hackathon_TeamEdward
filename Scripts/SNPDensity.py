from pybedtools import BedTool
import numpy as np
import matplotlib.pyplot as plt

with open('/mnt/disks/data-bam-tumor/Homo_sapiens_assembly19.fasta.fai','r') as f:
    chromosomes = {line.split('\t')[0]:int(line.split('\t')[1]) for line in f.readlines()[0:24]}

print chromosomes

chrom_Intervals = {k:np.append(np.arange(0,chromosomes[k],10000),chromosomes[k]) for k in chromosomes.keys()}

for k in chromosomes.keys():
    print chrom_Intervals[k][0:10]


with open('humanBedInterval.bed','w') as f:
    for k in chromosomes.keys():
        f.write('\n'.join('\n'.join('%s\t%d\t%d'%(k,chrom_Intervals[k][i],chrom_Intervals[k][i+1]-1) for i in range(len(chrom_Intervals[k])-1)) for k in chromosomes.keys())+'\n')

a = BedTool('humanBedInterval.bed')
b = BedTool('/mnt/disks/data-vcf/GSN79Tumor_normal.vcf')
print a.head()
print b.head()

a.coverage(b,hist=True).saveas('VCFCoverage.bed')

#for k in chromosomes.keys():
