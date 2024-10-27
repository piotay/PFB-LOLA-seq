# **LOLASEQ**

## **Developing pipelines for single cell RNA sequencing (scRNAseq) and, Assay for Transposase-Accessible Chromatin (ATAC-Seq)**.

- First item scRNAseq and ATACseq are next-generation sequencing technologies.
- Second item ATACseq is to determine how accessible the chromatin is for the transcription machinary.
  
## The datasets used

- First item NCBI GEO accession: GSE218392, 
- Second item Human cell Targeting Immune-Fibroblast Crosstalk in Myocardial Infarction and Cardiac Fibrosis

## Pipeline: scanpy
```python
import os  
import numpy as np  
import pandas as pd  
import scanpy as sc 
import anndata as ad
import muon as mu
```
 # For ATACseq, in addition,
 ### Import a module with ATAC-seq-related functions
 ```python
from muon import atac as ac

  mdata = mu.read("data/pbmc10k.h5mu")
  mdata
```

![image](https://github.com/user-attachments/assets/4edc7488-37db-455f-b69d-7b5d7519ffb1)

```python
mdata.var_names_make_unique()
mdata
```

![image](https://github.com/user-attachments/assets/c549248f-e0da-4357-8da9-e396cd81b633)

```python
from muon import atac as ac
atac = mdata.mod['atac']
sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=False, inplace=True)
mu.pl.histogram(atac, ['n_genes_by_counts', 'total_counts'], linewidth=0)
```





