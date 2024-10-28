#!usr/bin/env python3

import sys

inputfile=sys.argv[1]
outputfile=sys.argv[2]
depththreshold=int(sys.argv[3])
lengththreshold=int(sys.argv[4])
peaks=[]

#ensure both arguments are present, else exit
if len(sys.argv) != 5:
    print(f'Usage: {sys.argv[0]} <input file name> <output file name> <depth threshold for peaks> <length threshold>')
    exit(1)


#input file format looks like:
# 'Chromosome' \t nt coordinate startstreak \t nt coordinate endstreak \t ATAC read count depth
with open(inputfile,'r') as file: #open input file
    in_peak=False #are we in a peak or not? start false
    for line in file:
        line=line.rstrip()
        base=line.split('\t') #split into list: ['chromosome',coordinate,count]
        if int(base[3]) > depththreshold and in_peak == False: #if count is 20 or more and we are out of a peak:
            in_peak=True #change to in peak
            start_coord=int(base[1]) #remember start of peak coordinate
        if int(base[3]) <= depththreshold and in_peak == True: #if count is less than 20 and we are in a peak:
            in_peak = False #we are no longer in the peak
            end_coord=int(base[1]) #remember end of peak coordinate
            peaks.append([base[0],start_coord,end_coord]) #add the peak to list of peaks ['chr1',startcoord,endcoord]

#filter out short peaks (100 or less bases)
filtered_peaks=[]
for peak in peaks:
    if peak[2]-peak[1] > lengththreshold: #if peak is 100 or more nts:
        filtered_peaks.append(peak) #add to filtered peaks list
    else: #skip if not
        continue


#use this if we want to use std.out > 
for peak in filtered_peaks:
    print(f'{peak[0]}\t{peak[1]}\t{peak[2]}')

#use this if we want an explicit file output
with open(outputfile,'w') as output:
    for peak in filtered_peaks:
        output.write(f'{peak[0]}\t{peak[1]}\t{peak[2]}\n')
output.close()