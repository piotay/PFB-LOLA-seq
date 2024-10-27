# LOLA-SEQ
## _Fixing the freaking files_

While we were downloading our files for linked scRNAseq and scATACseq, we realized that the authors did not provide the necessary files to run the scATACseq part. We needed the following:

- ATAC_fragments.tsv (supplied)
- matrix? (supplied)
- fragment index (not provided)
- peak annotations (not provided)

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
