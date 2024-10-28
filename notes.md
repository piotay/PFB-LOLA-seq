# LOLA-SEQ
## _Fixing the freaking files_

While we were downloading our files for linked scRNAseq and scATACseq, we realized that the authors did not provide the necessary files to run the scATACseq part. We needed the following:

- ATAC_fragments.tsv (supplied)
- matrix? (supplied)
- fragment index (not provided)
- peak annotations (not provided)

## Obtaining the Peak Annontations
The fragments file is composed of all read fragments which have already been aligned to the human genome. It is structured with the following 5 fields:

- chromosome#
- start_position
- end_position
- UMI
- read_count

From this starting material, we needed to generate a file that had the following fields:
- Peak location
- Closest gene
- Distance from closest gene
- Location descriptor

## Step 1: Converting read fragments into coverage
We used __Bedtools Genome Coverage__ as the first step. This tool takes all read locations, and lines them up along the length of the human genome. 

![Alt text](https://bedtools.readthedocs.io/en/latest/_images/genomecov-glyph.png)

We slightly altered our two input files so that they would be compatible with this tool. 

From the original fragments file, we needed to remove had several lines of comment and reads that mapped onto chromosome scaffolds: We wrote the system output to a new file:
>'grep v # fragments.tsv > fragments_modified.tsv
>'grep -v -e GL -e KI fragments_modified.tsv > fragments_modified.tsv

chr1    10078   10308   GGAATCTTCCGTAAAC-1      1
chr1    10079   10261   GGCTAGACAACCTAAT-1      1
chr1    10079   10309   TATGGGCGTAAACAAG-1      1
chr1    10084   10317   GGTGTTGTCCGCCTCA-1      1

The other file needed contained the length of each human chromosome in nucleotides, which specified the length we were mapping onto.

chr1	249250621
chr2	243199373
chr3	198022430
chr4	191154276

There are several options for output, and we opted for the output which collapses our overlapping reads based on location. We also did not include zeros to decrease our file size.
>'bedtools genomecov =dz -i fragments_modified.tsv -g human_chrom_sizes.txt > genome_cov.txt

## Step 2: Calling peaks from the coverage


## Step 3: 

- Download .gff3 file from Gencode that annotates every gene in the human genome
- Convert to .bed file using gff2bed from BEDOPS package
  >`gff2bed < input.gff3 > output.bed`









## Peak Annotations:

### First need annotation of entire genome
- Download .gff3 file from Gencode that annotates every gene in the human genome
- Convert to .bed file using gff2bed from BEDOPS package
  >`gff2bed < input.gff3 > output.bed`
- filter to only genes
  >`awk '{if ($8 == "gene") print}' file.bed`
- filter to only protein-coding genes
  >`awk '/gene_type=protein_coding/ {print $0}' file.bed`
- Use python script `polish_bed.py` to update the ENSEMBL ID with the gene name
- filter to just location and gene
  >`% cut -f 1-4`
- this file is used to annotate peaks, should look like:
  >`chr1	65418	71585	OR4F5`
  > `chr1	450739	451678	OR4F29`

recieve file from nathan, do `bedtools closest -d -a nate_file -b gene_annotations`

`format_closest_overlap.py` on the output of bedtools closest, then run `define_distance.py`

output from this file is sent back to nate for scATAC pipeline
