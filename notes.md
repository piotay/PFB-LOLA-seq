# LOLA-SEQ
## _Fixing the freaking files_

While we were downloading our files for linked scRNAseq and scATACseq, we realized that the authors did not provide the necessary files to run the scATACseq part. We needed the following:

- ATAC_fragments.tsv (supplied)
- matrix? (supplied)
- fragment index (not provided)
- peak annotations (not provided)

## Peak Annotations:

- Download .gff3 file from Gencode that annotates every gene in the human genome
- Convert to .bed file using gff2bed from BEDOPS package
- <code>$ gff2bed < input.gff3 > output.bed
- filter to only genes
- <code> awk '{if ($8 == "gene") print}' file.bed
- filter to only protein-coding genes
- <code> awk '/gene_type=protein_coding/ {print $0}' file.bed
- Use script to update the ENSEMBL ID with the gene name
- Use <code>% cut -f 1-4
- to return only the gene location and name
- this file is used to annotate peaks

recieve file from nathan, do bedtools closest -d -a nate_file -b gene_annotations 

format_closest_overlap.py on the output of bedtools closest, then run define_distance.py

output from this file is sent back to nate for scATAC pipeline
