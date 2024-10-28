# snRNASeq Analysis

This dataset is derived from a patient with ischemic cardiomyopathy, or heart failure. The tissue is obtained during a heart transplant when a patient is receiving a donor heart. The tissue is processed to obtain information about transcriptomic signatures and chromatin accessibility on a single-cellular level.

NCBI GEO accession: GSE218392
Targeting Immune-Fibroblast Crosstalk in Myocardial Infarction and Cardiac Fibrosis
The fltered barcode matrix in .h5 file had two modalities, RNA and ATAC. RNA was used for one half of the analysis and ATAC was used for the other half.


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

Principal components are then used to compute the neighborhood graph for cells. The Uniform MAnifold Approximation and Projection (UMAP) represents the reduced dimensionality accounting for the variance in gene expression between groups of cells.

![Alttext](https://raw.githubusercontent.com/rdalipo1/PFB-LOLA-seq/refs/heads/main/pfb_lola_seq/scrna-seq/rna_atac_figs/UMAP.png)

Each group is a cluster with grouped by similar gene expression. The next step would be to annotate the cell clusters. To do this, we look at the top differentially expressed genes for each cluster.

![Alttext](https://raw.githubusercontent.com/rdalipo1/PFB-LOLA-seq/refs/heads/main/pfb_lola_seq/scrna-seq/rna_atac_figs/TopDEgenes.png)

This file is used moving forward in the ATACseq analysis.

# Chromatin accessibilty
The ATACseqdata is processed in a very similar pipeline to the snRNAseq data, generally with pre-processing, normalization, dimensionality reduction, and clustering being performed on sequences of open chromatin structures rather than transcription level.

## Preprocessing and quality control
For quality control, instead of working with gene counts we are working with peaks located at genes. We start with a file that contains a matrix containing the peak height of each gene in all collected nuclei. We need again filter using biologically relevant metrics.

![Alttext

ATAC-specific QC involves looking at nucleosome signal and transcription start site (TSS) enrichment. We expect chromatin accessibility enriched around TSSs.

## Principal component analysis, dimensionality reduction, and clustering
The data is again normalized, scaled, selected for highly variable features.

After principal component analysis and dimensionality reduction, another umap is compiled from clustering the chromatin accessibilty data.

