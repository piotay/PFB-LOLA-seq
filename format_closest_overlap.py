#!/usr/bin/env python3
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file) as ifh, open(output_file, 'w') as wfh:
    for line in ifh:
        temp = line.split('\t')
        line = temp[0]+'_'+temp[1]+'_'+temp[2]+'\t'+temp[6]+'\t'+temp[7]       
        wfh.write(line)
