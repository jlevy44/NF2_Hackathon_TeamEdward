import numpy as np
from pybedtools import BedTool
from collections import defaultdict
from cPickle import *
from random import randrange
import sys
# Test data


def genDataset(genes,testTrain): # second argument is test or train bed dictionary
    dataset = {'SNP':defaultdict(list),'indel':defaultdict(list)}
    #random_index = randrange(0,len(genes))
    #with open('out.txt','w') as f:
    if len(genes) > 1:
        for gene in genes[0:10]:#random_index[0:11]:
            #gene = genes[i]
            print gene
            if gene and gene.startswith('1\t') or gene.startswith('22\t'):
                geneInfo = gene.split('\t')
                interval = map(int,geneInfo[1:3])
                #f.write('\t'.join([geneInfo[0]]+geneInfo[1:3])+'\n')
                #bin1 = np.arange(interval[0],interval[1],100)
                #bin2 = np.arange(interval[0]+50,interval[1],100)
                #geneBed = BedTool(gene,from_string=True)
                #SNPtest = testTrain['SNP'].intersect(geneBed,wa=True)
                #indelTest = testTrain['indel'].intersect(geneBed,wa=True)
                #snpDensity = []
                #indelDensity = []
                #print [np.arange(interval[0],interval[1],100),np.arange(interval[0]+50,interval[1],100)]
                for bin in [np.arange(interval[0],interval[1],100),np.arange(interval[0]+50,interval[1],100)]:
                    for i in range(len(bin)-1):
                        interval = bin[i:i+2]
                        geneNaming = geneInfo[3].strip('\n')+'|'+'-'.join([geneInfo[0]]+map(str,interval))#geneInfo[0:3]
                        #f.write('Gene Name: ' + geneNaming + '\n')
                        #try:
                        densityBedInt = BedTool('\n'.join(np.vectorize(lambda x: geneInfo[0]+'\t%d\t%d'%(x-10,x+10))(np.arange(interval[0]+10,interval[1]-10,5))),from_string=True)
                        #except:
                        #    print '\n'.join(np.vectorize(lambda x: geneInfo[0]+'\t%d\t%d'%(x-5,x+5))(np.arange(interval[0]+5,interval[1]-5)))

                        #try:
                        densitySNP = np.vectorize(lambda line: float(line.split('\t')[-1]))(filter(None,str(densityBedInt.coverage(testTrain['SNP'])).split('\n')))
                        #except:
                        #    print str(densityBedInt.coverage(testTrain['SNP'])).split('\n')
                        densityIndel = np.vectorize(lambda line: float(line.split('\t')[-1]))(filter(None,str(densityBedInt.coverage(testTrain['indel'])).split('\n')))
                        dataset['SNP'][geneNaming] = densitySNP
                        dataset['indel'][geneNaming] = densityIndel
                        #f.write(geneNaming+'\n')
        #f.write('FINISH 1\n')#testTrain['SNP'].head()
    else:
        for gene in genes:
            print gene
            if gene and gene.startswith('1\t') or gene.startswith('22\t'):
                geneInfo = gene.split('\t')
                interval = map(int,geneInfo[1:3])
                #f.write('\t'.join([geneInfo[0]]+geneInfo[1:3])+'\n')
                #bin1 = np.arange(interval[0],interval[1],100)
                #bin2 = np.arange(interval[0]+50,interval[1],100)
                #geneBed = BedTool(gene,from_string=True)
                #SNPtest = testTrain['SNP'].intersect(geneBed,wa=True)
                #indelTest = testTrain['indel'].intersect(geneBed,wa=True)
                #snpDensity = []
                #indelDensity = []
                #print [np.arange(interval[0],interval[1],100),np.arange(interval[0]+50,interval[1],100)]
                for bin in [np.arange(interval[0],interval[1],100),np.arange(interval[0]+50,interval[1],100)]:
                    for i in range(len(bin)-1):
                        interval = bin[i:i+2]
                        geneNaming = geneInfo[3].strip('\n')+'|'+'-'.join([geneInfo[0]]+map(str,interval))#geneInfo[0:3]
                        #f.write('Gene Name: ' + geneNaming + '\n')
                        #try:
                        densityBedInt = BedTool('\n'.join(np.vectorize(lambda x: geneInfo[0]+'\t%d\t%d'%(x-10,x+10))(np.arange(interval[0]+10,interval[1]-10,5))),from_string=True)
                        #except:
                        #    print '\n'.join(np.vectorize(lambda x: geneInfo[0]+'\t%d\t%d'%(x-5,x+5))(np.arange(interval[0]+5,interval[1]-5)))

                        #try:
                        densitySNP = np.vectorize(lambda line: float(line.split('\t')[-1]))(filter(None,str(densityBedInt.coverage(testTrain['SNP'])).split('\n')))
                        #except:
                        #    print str(densityBedInt.coverage(testTrain['SNP'])).split('\n')
                        densityIndel = np.vectorize(lambda line: float(line.split('\t')[-1]))(filter(None,str(densityBedInt.coverage(testTrain['indel'])).split('\n')))
                        dataset['SNP'][geneNaming] = densitySNP
                        dataset['indel'][geneNaming] = densityIndel
        print 'FINISH 1'
        return dataset






trainBed = {'SNP':BedTool('train_SNP.bed'),'indel':BedTool('train_indel.bed')}
testBed = {'SNP':BedTool('test_SNP.bed'),'indel':BedTool('test_indel.bed')}

BedTool('human_genes.bed').sort().intersect(testBed['SNP'].cat(testBed['indel']),wa=True).sort().merge(c=4,o='distinct').saveas('human_genes_Check.bed')
trainGenes = BedTool('human_genes.bed').sort().intersect(trainBed['SNP'].cat(trainBed['indel']),wa=True).sort().merge(c=4,o='distinct').saveas('trainGenes.bed')

with open('trainGenes.bed','r') as f:
    genes = f.readlines()
#generate training data set

trainData = genDataset(genes,trainBed)

dump(trainData,open('trainData.p','wb'))



# generate final data set
humanBed = 'human_genes_Check.bed'
with open(humanBed,'r') as f:
    genes = f.readlines()

dump(genDataset(genes,testBed),open('testData.p','wb'))

