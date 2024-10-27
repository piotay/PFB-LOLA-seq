# **LOLASEQ**

## **Developing pipelines for single cell RNA sequencing (scRNAseq) and, Assay for Transposase-Accessible Chromatin (ATAC-Seq)**.

- First item scRNAseq and ATACseq are next-generation sequencing technologies.
- Second item ATACseq is to determine how accessible the chromatin is for the transcription machinary.
  
## The datasets used

- First item NCBI GEO accession: GSE218392, 
- Second item Human cell Targeting Immune-Fibroblast Crosstalk in Myocardial Infarction and Cardiac Fibrosis

## Pipeline: scanpy

`import os`  
`import numpy as np`  
`import pandas as pd`  
`import scanpy as sc`  
`import anndata as ad`  
`import muon as mu`

 #For ATACseq, in addition,
 ### Import a module with ATAC-seq-related functions
 from muon import atac as ac

mdata = mu.read("data/pbmc10k.h5mu")
mdata






