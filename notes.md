# LOLA-SEQ
## Fixing the freaking files

While we were downloading our files for linked scRNAseq and scATACseq, we realized that the authors did not provide the necessary files to run the scATACseq part. We needed the following:

- ATAC_fragments.tsv (supplied)
- Feature_barcode_matrix.h5 (supplied)
- fragment_index.tsv (not provided)
- peak_annotations.tsv (not provided)

## Obtaining the peak annotations
The fragments file is composed of all read fragments which have already been aligned to the human genome. It is structured with the following 5 fields:
- Chromosome #
- Start nt position
- End nt position
- UMI
- Read

From this starting material, we needed to generate a peak annotations file that had the following fields:
- Peak location on chromosome
- Closest gene
- Distance from closest gene
- Location descriptor

![Alt text](https://www.basepairtech.com/wp-content/uploads/2022/06/ATAC_USP47.png)

## Step 1: Converting read fragments into coverage
We used __Bedtools Genome Coverage__ as the first step. This tool takes all read locations, and lines them up along the length of the human genome. 

![Alt text](https://bedtools.readthedocs.io/en/latest/_images/genomecov-glyph.png)

We slightly altered our two input files so that they would be compatible with this tool. 

From the original fragments file, we needed to remove had several lines of comment and reads that mapped onto chromosome scaffolds: We wrote the system output to a new file:
>'grep v # fragments.tsv > fragments_modified.tsv
>'grep -v -e GL -e KI fragments_modified.tsv > fragments_modified.tsv

 ```
  chr1  10078  10308  GGAATCTTCCGTAAAC-1  1 
  chr1  10079  10261  GGCTAGACAACCTAAT-1  1 
  chr1  10079  10309  TATGGGCGTAAACAAG-1  1 
  chr1  10084  10317  GGTGTTGTCCGCCTCA-1  1 
```
The other file needed contained the length of each human chromosome in nucleotides, which specified the length we were mapping onto.
```
chr1	249250621
chr2	243199373
chr3	198022430
chr4	191154276
```
There are several options for output, and we opted for the output which collapses our overlapping reads based on location. We also did not include zeros to decrease our file size.
>'bedtools genomecov =dz -i fragments_modified.tsv -g human_chrom_sizes.txt > genome_cov.txt`

## Step 2: Calling peaks from the coverage


## Step 3: Generating Peak Annotations

### Step 3.1: Finding human genome annotations

- We first downloaded a .gff3 file from Gencode that annotates every gene in the human genome
- Next, we converted this .gff3 file to a .bed file using gff2bed from BEDOPS package
  ```
  $ gff2bed < input.gff3 > output.bed
  ```
- filter to only genes
  ```
  $ awk '{if ($8 == "gene") print}' file.bed
  ```
- filter to only protein-coding genes
  ```
  $ awk '/gene_type=protein_coding/ {print $0}' file.bed
  ```
- Use python script `polish_bed.py` to update the ENSEMBL ID with the gene name
- filter to just location and gene
  ```
  $ cut -f 1-4
  ```
- this file is used to annotate peaks, should look like:
  ```
  chr1	65418	71585	OR4F5
  chr1	450739	451678	OR4F29
  ```

### Step 3.2: Annotate the Peaks
Using the peaks generated from step 2, we determined which genes these peaks are closest to using BedTools
```
$ bedtools closest -d -a <peak_calls file> -b <gene_annotations_file>
```

We next wrote two python scripts to help format the output from the `bedtoosl closest` tool.
We used `format_closest_overlap.py` on the output of bedtools closest, then `define_distance.py` to call the peaks as either "distal" from the gene or "promoter" if the peak is directly overlapping the called gene.

The output from this file is the final peak annotations file, and is the final piece of our input for the scATAC pipeline!
