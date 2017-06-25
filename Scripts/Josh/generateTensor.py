import numpy as np
from pybedtools import BedTool
from collections import defaultdict
from cPickle import *
import sys
# Test data


def genDataset(genes,testTrain): # second argument is test or train bed dictionary
    dataset = {'SNP':defaultdict(list),'indel':defaultdict(list)}
    with open('Error.txt','w') as f:
        for gene in genes:
            if gene:
                try:
                    geneInfo = gene.split('\t')
                    interval = map(int,geneInfo[1:3])
                    #bin1 = np.arange(interval[0],interval[1],100)
                    #bin2 = np.arange(interval[0]+50,interval[1],100)
                    #geneBed = BedTool(gene,from_string=True)
                    #SNPtest = testTrain['SNP'].intersect(geneBed,wa=True)
                    #indelTest = testTrain['indel'].intersect(geneBed,wa=True)
                    #snpDensity = []
                    #indelDensity = []
                    for bin in [np.arange(interval[0],interval[1],100),np.arange(interval[0]+50,interval[1],100)]:
                        for i in range(len(bin)):
                            interval = bin[i:i+2]
                            geneNaming = geneInfo[3]+'|'+'-'.join(map(str,interval))#geneInfo[0:3]
                            densityBedInt = BedTool('\n'.join(np.vectorize(lambda x: geneInfo[0]+'\t%d\t%d'%(x-5,x+5))(np.arange(interval[0]+5,interval[1]-5))),from_string=True)
                            densitySNP = np.vectorize(lambda line: int(line.split('\t')[-1]))(str(densityBedInt.coverage(testTrain['SNP'])).split('\n'))
                            densityIndel = np.vectorize(lambda line: int(line.split('\t')[-1]))(str(densityBedInt.coverage(testTrain['indel'])).split('\n'))
                            dataset['SNP'][geneNaming] = densitySNP
                            dataset['indel'][geneNaming] = densityIndel
                except:
                    f.write(gene+' Error:'+str(sys.exc_info()[0])+'\n')
    return dataset






trainBed = {'SNP':BedTool('test_SNP.bed'),'indel':BedTool('test_indel.bed')}
testBed = {'SNP':BedTool('train_SNP.bed'),'indel':BedTool('train_indel.bed')}


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

