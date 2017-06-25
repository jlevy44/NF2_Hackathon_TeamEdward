with open('human.bed','r') as f:
    lines = f.readlines()

with open('human_genes.bed','w') as f:
    f.writelines(['\t'.join(line.split('\t')[0:4]+[line.split('\t')[9].split(';')[0].replace('ID=gene:','')]) for line in lines if line and line.split('\t')[7]=='gene'])