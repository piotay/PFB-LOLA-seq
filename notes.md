# LOLA-SEQ
# Goals
This dataset is derived from a patient with ischemic cardiomyopathy, or heart failure. The tissue is obtained during a heart transplant when a patient is receiving a donor heart. The tissue is processed to obtain information about transcriptomic signatures and chromatin accessibility on a single-cellular level.

__NCBI GEO accession: GSE218392__

__Targeting Immune-Fibroblast Crosstalk in Myocardial Infarction and Cardiac Fibrosis__

The filtered barcode matrix in .h5 file had two modalities, RNA and ATAC. RNA was used for one half of the analysis and ATAC was used for the other half.

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
```
$ grep v # fragments.tsv > fragments_modified.tsv
```
```
$ grep -v -e GL -e KI fragments_modified.tsv > fragments_modified.tsv
```
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
```
$ bedtools genomecov =dz -i fragments_modified.tsv -g human_chrom_sizes.txt > genome_cov.txt
```

## Step 2: Calling peaks from the coverage


## Step 3: Generating Peak Annotations

### Step 3.1: Finding human genome annotations

- We first downloaded a .gff3 file from Gencode that annotates every gene in the human genome
  ```
  $ head gencode.v47.gff3
  ```
  ```
  ##gff-version 3
  #description: evidence-based annotation of the human genome (GRCh38), version 47 (Ensembl 113)
  #provider: GENCODE
  #contact: gencode-help@ebi.ac.uk
  #format: gff3
  #date: 2024-07-19
  ##sequence-region chr1 1 248956422
  chr1	HAVANA	gene	11121	24894	.	+	.	ID=ENSG00000290825.2;gene_id=ENSG00000290825.2;gene_type=lncRNA;gene_name=DDX11L16;level=2;tag=overlaps_pseudogene
  chr1	HAVANA	transcript	11426	14409	.	+	.	ID=ENST00000832828.1;Parent=ENSG00000290825.2;gene_id=ENSG00000290825.2;transcript_id=ENST00000832828.1;gene_type=lncRNA;gene_name=DDX11L16;transcript_type=lncRNA;transcript_name=DDX11L16-264;level=2;tag=basic,TAGENE
  chr1	HAVANA	exon	11426	11671	.	+	.	ID=exon:ENST00000832828.1:1;Parent=ENST00000832828.1;gene_id=ENSG00000290825.2;transcript_id=ENST00000832828.1;gene_type=lncRNA;gene_name=DDX11L16;transcript_type=lncRNA;transcript_name=DDX11L16-264;exon_number=1;exon_id=ENSE00004248702.1;level=2;tag=basic,TAGENE
  chr1	HAVANA	exon	12010	12227	.	+	.	ID=exon:ENST00000832828.1:2;Parent=ENST00000832828.1;gene_id=ENSG00000290825.2;transcript_id=ENST00000832828.1;gene_type=lncRNA;gene_name=DDX11L16;transcript_type=lncRNA;transcript_name=DDX11L16-264;exon_number=2;exon_id=ENSE00004248735.1;level=2;tag=basic,TAGENE
  chr1	HAVANA	exon	12613	12721	.	+	.	ID=exon:ENST00000832828.1:3;Parent=ENST00000832828.1;gene_id=ENSG00000290825.2;transcript_id=ENST00000832828.1;gene_type=lncRNA;gene_name=DDX11L16;transcript_type=lncRNA;transcript_name=DDX11L16-264;exon_number=3;exon_id=ENSE00003582793.1;level=2;tag=basic,TAGENE
  chr1	HAVANA	exon	13221	14409	.	+	.	ID=exon:ENST00000832828.1:4;Parent=ENST00000832828.1;gene_id=ENSG00000290825.2;transcript_id=ENST00000832828.1;gene_type=lncRNA;gene_name=DDX11L16;transcript_type=lncRNA;transcript_name=DDX11L16-264;exon_number=4;exon_id=ENSE00004248703.1;level=2;tag=basic,TAGENE
  ```
- Next, we converted this .gff3 file to a .bed file using gff2bed from BEDOPS package
  ```
  $ gff2bed < input.gff3 > output.bed
  ```
- Since the .gff3 file contains exons as well as genes, we filtered to only include genes
  ```
  $ awk '{if ($8 == "gene") print}' file.bed
  ```
- Further, we filtered to only protein-coding genes
  ```
  $ awk '/gene_type=protein_coding/ {print $0}' file.bed
  ```
  ```
  chr1	65418	71585	ENSG00000186092.7	.	+	HAVANA	gene	.ID=ENSG00000186092.7;gene_id=ENSG00000186092.7;gene_type=protein_coding;gene_name=OR4F5;level=2;hgnc_id=HGNC:14825;havana_gene=OTTHUMG00000001094.4
  chr1	450739	451678	ENSG00000284733.2	.	-	HAVANA	gene	.ID=ENSG00000284733.2;gene_id=ENSG00000284733.2;gene_type=protein_coding;gene_name=OR4F29;level=2;hgnc_id=HGNC:31275;havana_gene=OTTHUMG00000002860.3
  chr1	685715	686654	ENSG00000284662.2	.	-	HAVANA	gene	.ID=ENSG00000284662.2;gene_id=ENSG00000284662.2;gene_type=protein_coding;gene_name=OR4F16;level=2;hgnc_id=HGNC:15079;havana_gene=OTTHUMG00000002581.3
  ```
- Use python script `polish_bed.py` to update the ENSEMBL ID with the gene name
  ```
  chr1	65418	71585	OR4F5	.	+	HAVANA	gene	.	ID=ENSG00000186092.7;gene_id=ENSG00000186092.7;gene_type=protein_coding;gene_name=OR4F5;level=2;hgnc_id=HGNC:14825;havana_gene=OTTHUMG00000001094.4
  chr1	450739	451678	OR4F29	.	-	HAVANA	gene	.	ID=ENSG00000284733.2;gene_id=ENSG00000284733.2;gene_type=protein_coding;gene_name=OR4F29;level=2;hgnc_id=HGNC:31275;havana_gene=OTTHUMG00000002860.3
  chr1	685715	686654	OR4F16	.	-	HAVANA	gene	.	ID=ENSG00000284662.2;gene_id=ENSG00000284662.2;gene_type=protein_coding;gene_name=OR4F16;level=2;hgnc_id=HGNC:15079;havana_gene=OTTHUMG00000002581.3
  ```
- filter to just location and gene
  ```
  $ cut -f 1-4
  ```
  ```
  chr1	65418	71585	OR4F5
  chr1	450739	451678	OR4F29
  chr1	685715	686654	OR4F16
  ```
- This is the file we used to annotate peaks!

### Step 3.2: Annotate the Peaks
Using the peaks generated from step 2, we determined which genes these peaks are closest to using BedTools
```
$ bedtools closest -d -a <peak_calls file> -b <gene_annotations_file>
```
The -d option will include the distance from the beak to the annotated gene.

This generated the following file:
```
chr1	10103	10339	chr1	65418	71585	OR4F5	55080
chr1	181061	181163	chr1	65418	71585	OR4F5	109477
chr1	181357	181627	chr1	65418	71585	OR4F5	109773
chr1	191297	191679	chr1	65418	71585	OR4F5	119713
chr1	631501	631696	chr1	685715	686654	OR4F16	54020
chr1	634593	634715	chr1	685715	686654	OR4F16	51001
chr1	777624	778145	chr1	685715	686654	OR4F16	90971
```

We next wrote two python scripts to help format the output from the `bedtoosl closest` tool.
We used `format_closest_overlap.py` on the output of bedtools closest to generate the following output:
```
chr1_10103_10339	OR4F5	55080
chr1_181061_181163	OR4F5	109477
chr1_181357_181627	OR4F5	109773
chr1_191297_191679	OR4F5	119713
chr1_631501_631696	OR4F16	54020
chr1_634593_634715	OR4F16	51001
chr1_777624_778145	OR4F16	90971
```
Then used `define_distance.py` to call the peaks as either "distal" from the gene or "promoter" if the peak is directly overlapping the called gene.
The output from this file (seen below) is the final peak annotations file, and is the final piece of our input for the scATAC pipeline!
```
chr1_10103_10339	OR4F5	55080	distal
chr1_181061_181163	OR4F5	109477	distal
chr1_181357_181627	OR4F5	109773	distal
chr1_191297_191679	OR4F5	119713	distal
chr1_631501_631696	OR4F16	54020	distal
chr1_634593_634715	OR4F16	51001	distal
chr1_777624_778145	OR4F16	90971	distal
```



# snRNAseq
![Alttext](https://cdn.10xgenomics.com/image/upload/v1709930681/blog/GEM-X%20Launch%20blog/Figure_1.png)

__Libraries used:__
- scanpy: a Pyhton toolkit for analyzing single-cell gene expression data
- muon: a Python framework designed to work with multimodal omics data
- pandas
- numpy
- matlabplot

## Preprocessing and quality control
We start with a file that contains a matrix containing the level of expression of each gene in all collected nuclei. We need to filter using biologically relevant metrics.

![Alttext](https://raw.githubusercontent.com/rdalipo1/PFB-LOLA-seq/refs/heads/main/pfb_lola_seq/scrna-seq/rna_atac_figs/QC_prefilter.png)


```
#filter by minimum and maximum gene counts
mu.pp.filter_obs(rna, 'n_genes_by_counts', lambda x: (x >= 200) & (x < 4000))
#filter by percentage mitochondria to remove potential dead or stressed cells
mu.pp.filter_obs(rna, 'pct_counts_mt', lambda x: x < 25)
rna = rna[rna.obs['percent_ribo'] < 0.05, :]
```

![Alttext](https://raw.githubusercontent.com/rdalipo1/PFB-LOLA-seq/refs/heads/main/pfb_lola_seq/scrna-seq/rna_atac_figs/QC_postfilter.png)

Genes and gene counts between cells are normalized, scaled, and selected for variability. 

## Dimensionality reduction and clustering
Principal component analysis returns principal compenents that describe some measure of the varibaility in the data. 

![Alttext](https://raw.githubusercontent.com/rdalipo1/PFB-LOLA-seq/refs/heads/main/pfb_lola_seq/scrna-seq/rna_atac_figs/PCA_elbowplot.png)

Principal components are then used to compute the neighborhood graph for cells. The Uniform Manifold Approximation and Projection (UMAP) represents the reduced dimensionality accounting for the variance in gene expression between groups of cells.

![Alttext](https://raw.githubusercontent.com/rdalipo1/PFB-LOLA-seq/refs/heads/main/pfb_lola_seq/scrna-seq/rna_atac_figs/UMAP.png)

Each group is a cluster with grouped by similar gene expression. The next step would be to annotate the cell clusters. To do this, we look at the top differentially expressed genes for each cluster.

![Alttext](https://raw.githubusercontent.com/rdalipo1/PFB-LOLA-seq/refs/heads/main/pfb_lola_seq/scrna-seq/rna_atac_figs/TopDEgenes.png)

This file is used moving forward in the ATACseq analysis.

# Chromatin accessibilty
The ATACseqdata is processed in a very similar pipeline to the snRNAseq data, generally with pre-processing, normalization, dimensionality reduction, and clustering being performed on sequences of open chromatin structures rather than transcription level.

## Preprocessing and quality control
For quality control, instead of working with gene counts we are working with peaks located at genes. We start with a file that contains a matrix containing the peak height of each gene in all collected nuclei. We need again filter using biologically relevant metrics.

![Alttext](https://raw.githubusercontent.com/piotay/PFB-LOLA-seq/refs/heads/main/pfb_lola_seq/scrna-seq/atac/preprocessing.png)

![Alttext](https://raw.githubusercontent.com/piotay/PFB-LOLA-seq/refs/heads/main/pfb_lola_seq/scrna-seq/atac/postprocessing.png)

ATAC-specific QC involves looking at nucleosome signal and transcription start site (TSS) enrichment. The nucleosome signal is the ratio of mono-nucleosomal to nucleosome-free fragments, also used as a signal-to-noise ratio in each cell. We expect chromatin accessibility enriched around TSSs.

![Alttext](https://raw.githubusercontent.com/piotay/PFB-LOLA-seq/refs/heads/main/pfb_lola_seq/scrna-seq/atac/nucleosomesignal.png)

![Alttext](https://raw.githubusercontent.com/piotay/PFB-LOLA-seq/refs/heads/main/pfb_lola_seq/scrna-seq/atac/TSSenrichment.png)

## Principal component analysis, dimensionality reduction, and clustering
The data is again normalized, scaled, selected for highly variable features.

After principal component analysis and dimensionality reduction, another umap is compiled from clustering the chromatin accessibilty data.

![Alttext](https://raw.githubusercontent.com/piotay/PFB-LOLA-seq/refs/heads/main/pfb_lola_seq/scrna-seq/atac/UMAP.png)
![Alttext](https://raw.githubusercontent.com/piotay/PFB-LOLA-seq/refs/heads/main/pfb_lola_seq/scrna-seq/atac/dotplot.png)

