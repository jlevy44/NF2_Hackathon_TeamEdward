from pybedtools import BedTool
import numpy as np

with open('/mnt/disks/data-bam-tumor/Homo_sapiens_assembly19.fasta.fai','r') as f:
    chromosomes = {line.split('\t')[0]:int(line.split('\t')[1]) for line in f.readlines()[0:25]}

chrom_Intervals = {k:np.append(np.arange(0,chromosomes[k],10000),chromosomes[k]) for k in chromosomes.keys()}


with open('humanBedInterval.bed','w') as f:
    for k in chromosomes.keys():
        f.write('\n'.join('%s\t%d\t%d'%(k,chrom_Intervals[i],chrom_Intervals[i+1]-1) for i in range(len(chrom_Intervals)-1))+'\n')

a = BedTool('humanBedInterval.bed')
b = BedTool('/mnt/disks/data-vcf/GSN79Tumor_normal.vcf')
print a.head()
print b.head()

a.coverage(b,hist=True).saveas('VCFCoverage.bed')