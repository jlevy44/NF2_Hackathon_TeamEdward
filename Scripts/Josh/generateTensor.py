import numpy as np
from pybedtools import BedTool
from collections import defaultdict
from cPickle import *
# Test data


def genDataset(genes,testTrain): # second argument is test or train bed dictionary
    dataset = defaultdict(list)
    for gene in genes:
        if gene:
            try:
                geneInfo = gene.split('\t')
                geneNaming = geneInfo[3]+'|'+'-'.join(geneInfo[0:3])
                interval = map(int,geneInfo[1:3])
                #geneBed = BedTool(gene,from_string=True)
                #SNPtest = testTrain['SNP'].intersect(geneBed,wa=True)
                #indelTest = testTrain['indel'].intersect(geneBed,wa=True)
                #snpDensity = []
                #indelDensity = []
                densityBedInt = BedTool('\n'.join(np.vectorize(lambda x: geneInfo[0]+'\t%d\t%d'%(x-5,x+5))(np.arange(interval[0]+5,interval[1]-5))),from_string=True)
                densitySNP = np.vectorize(lambda line: int(line.split('\t')[-1]))(str(densityBedInt.coverage(testTrain['SNP'])).split('\n'))
                densityIndel = np.vectorize(lambda line: int(line.split('\t')[-1]))(str(densityBedInt.coverage(testTrain['indel'])).split('\n'))
                dataset[geneNaming] = [densitySNP,densityIndel]
            except:
                print gene
    return dataset






trainBed = {'SNP':BedTool('file'),'indel':BedTool('file')}
testBed = {'SNP':BedTool('file'),'indel':BedTool('file')}


trainGenes = BedTool('human_genes.bed').sort().intersect(trainBed['SNP'].cat(trainBed['indel']),wa=True).saveas('trainGenes.bed')

with open('trainGenes.bed','r') as f:
    genes = f.readlines()
#generate training data set
dump(genDataset(genes,trainBed),open('trainData.p','wb'))



# generate final data set
humanBed = 'human_genes.bed'
with open(humanBed,'r') as f:
    genes = f.readlines()

dump(genDataset(genes,testBed),open('testData.p','wb'))

