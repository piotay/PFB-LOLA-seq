#!/usr/bin/env/python3
import re

'''
Script that takes in an input file (gencode human genome, filtered to only protein coding genes) and changes column 3 (ENSEMBL ID) to the gene name based on same annotation

Outputs fixed file, not filtered any further
'''



ifh = 'protein_coding_gencode.bed'
ofh = 'gene_protein_coding_gencode.bed'

with open(ifh) as file, open(ofh, 'w') as write_file:
    for line in file:
        line = line.rstrip()
        result = re.search(r'gene_name=(.*?);', line)
        line = re.sub(r'\sENSG\S+', '\t'+result.group(1), line)
        write_file.write(line+'\n')
        

