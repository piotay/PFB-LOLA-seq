#!/usr/bin/env python3
'''
Function that takes a .bed file as an input, formatted as such:
chr1_19494_29384 \t gene \t distance 

outputs file that adds a column calling the peak (based on distance)  as either: 

overlaps 1 gene: promoter
overlaps 0 genes: distal

'''
import pandas as pd
import numpy as np
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

df = pd.read_csv(input_file, sep='\t', header=None)
df['description'] = np.where(df[2] > 0, 'distal', 'promoter')

df.to_csv(output_file, sep='\t', header=False, index=False)
