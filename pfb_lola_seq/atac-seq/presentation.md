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
![image](https://github.com/user-attachments/assets/d58edd9c-4f0a-41e9-bc20-969d2ff43824)

```python
mu.pp.filter_var(atac, 'n_cells_by_counts', lambda x: x >= 10)
print(f"Before: {atac.n_obs} cells")
mu.pp.filter_obs(atac, 'total_counts', lambda x: (x >= 1000) & (x <= 80000))
print(f"(After total_counts: {atac.n_obs} cells)")
mu.pp.filter_obs(atac, 'n_genes_by_counts', lambda x: (x >= 100) & (x <= 30000))
print(f"After: {atac.n_obs} cells")
```
<span style="font-size: larger;">Before: 11909 cells</span>  
<span style="font-size: larger;">(After total_counts: 11564 cells)</span>  
<span style="font-size: larger;">After: 11564 cells</span>

```python
mu.pl.histogram(atac, ['n_genes_by_counts', 'total_counts'], linewidth=0)
```
![image](https://github.com/user-attachments/assets/2038603f-942f-455b-a398-40aeb83058f8)

```python
ac.pl.fragment_histogram(atac, region='chr1:1-2000000')
```
![image](https://github.com/user-attachments/assets/7fb78d6f-efdc-4ab8-a958-7b48faef0b9d)

```python
ac.tl.nucleosome_signal(atac, n=1e6)
mu.pl.histogram(atac, "nucleosome_signal", linewidth=0)
```
![image](https://github.com/user-attachments/assets/909a176a-6086-4c5d-833d-7c20c50e050a)

```python
ac.tl.get_gene_annotation_from_rna(mdata['rna']).head(3)
```














