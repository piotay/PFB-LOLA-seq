# **LOLASEQ**

## **Developing pipelines for single cell RNA sequencing (scRNAseq) and, Assay for Transposase-Accessible Chromatin (ATAC-Seq)**.

 **scRNAseq** and **ATACseq** are next-generation sequencing technologies.
- **ATACseq** is used to determine how accessible the chromatin is for the transcription machinery.
  
## The datasets used

**NCBI GEO accession:** GSE218392  
**Human cell Targeting Immune-Fibroblast Crosstalk in Myocardial Infarction and Cardiac Fibrosis**

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
| Chromosome | Start | End   | gene_id         | gene_name  |
|------------|-------|-------|-----------------|------------|
| MIR1302-2HG | chr1  | 29553 | 30267 | ENSG00000243485 | MIR1302-2HG |
| FAM138A     | chr1  | 36080 | 36081 | ENSG00000237613 | FAM138A     |
| OR4F5       | chr1  | 65418 | 69055 | ENSG00000186092 | OR4F5       |


```python
tss = ac.tl.tss_enrichment(mdata, n_tss=1000)
```
![image](https://github.com/user-attachments/assets/723e0148-d6f7-441d-ab71-27d72db6a559)

```python
tss
```
![image](https://github.com/user-attachments/assets/37ca0304-2c94-48c4-a103-9f184c73c513)

```python
ac.pl.tss_enrichment(tss)
```

![image](https://github.com/user-attachments/assets/38338970-8aa9-4eae-824c-0fbfa9222e62)

```python
atac.layers["counts"] = atac.X.copy()
sc.pp.normalize_total(atac, target_sum=1e4)
sc.pp.log1p(atac)
atac.layers["lognorm"] = atac.X.copy()
sc.pp.highly_variable_genes(atac, min_mean=0.05, max_mean=1.5, min_disp=.5)
sc.pl.highly_variable_genes(atac)
```

![image](https://github.com/user-attachments/assets/3db30301-ca8a-4cc5-97b2-fdba2ed43cab)

```python
np.sum(atac.var.highly_variable)
```

'np.int64(11771)'

```python
sc.pp.scale(atac, max_value=10)
sc.tl.pca(atac, svd_solver='arpack')
ac.pl.pca(atac, color=['NRCAM', 'SLC1A2', 'SRGN', 'VCAN'], layer='lognorm', func='mean')
```

![image](https://github.com/user-attachments/assets/036561ae-d4b5-4da3-befe-98437827b959)

```python
sc.pl.pca_variance_ratio(atac, log=True)
```

![image](https://github.com/user-attachments/assets/e7f99d4e-e330-4ac9-aedd-a3eb9e58213f)

```python
sc.pp.neighbors(atac, n_neighbors=10, n_pcs=20)
sc.tl.leiden(atac, resolution=.5)
```
![image](https://github.com/user-attachments/assets/810d965d-2fff-45a0-bdde-4e67eb6bdfd5)


```python
sc.tl.umap(atac, spread=1., min_dist=.5, random_state=11)
sc.pl.umap(atac, color="leiden", legend_loc="on data")
```

![image](https://github.com/user-attachments/assets/4568c0e4-6d88-4f38-bf0b-ea33a338b104)

```python
ac.tl.rank_peaks_groups(atac, 'leiden', method='t-test')
```

![image](https://github.com/user-attachments/assets/0a981e17-54db-4fc0-8210-c46994ea4432)

```python
result = atac.uns['rank_genes_groups']
groups = result['names'].dtype.names

pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'genes', 'pvals']}).head(10)
```
![image](https://github.com/user-attachments/assets/87eecff2-0109-408c-b59a-3dddfe0076c6)
























