---
layout: default
title: Non-Nueroepithelial Data 
parent: Classifying single cell or nuclei RNA-seq
nav_order: 2 
---

# Input For Classification
The ccAFv2 classifier requires a thoroughly quality-controlled Seurat
object. SCTransform is preferred, standard normalization only applies to the
highly variable genes. This can exclude genes needed for the accurate
classification of the cell cycle. The PredictCellCycle
function re-runs SCTransform to retain all genes in the dataset. 

An example QC pipeline used on the U5 human neural stem cell (hNSC) dataset is available
[here](https://github.com/plaisier-lab/FIXLINK).

## Non-Neuroepithelial Input Data

The Human Skeletal Muscle Satellite Cells (HSkMSC) dataset used to test the ccAFv2 is
available
[here](https://zenodo.org/records/24207046/files/HSkMSC_cr_normalized_ensembl.rds?download=).
This data has been QC'd and normalized using SCTransform following our
best practices described above.
### Load packages 
```r
library(ccAFv2)
library(Seurat)
```

### Load in sample data
```r
seurat_obj = readRDS('hSkMSC_normalized_ensembl.rds')
```

## Predicting cell-cycle states
Examining the number of cell-cycle marker genes is crucial for accuracy. There are 861 cell cycle marker genes in ccAFv2 Predictions significantly worsen if fewer than 689 marker genes (80%) are present. Setting do_sctransform = TRUE maximizes marker gene overlap. Timing and 93/93 values may vary by dataset, which is normal.

```r
seurat_obj = PredictCellCycle(seurat_obj,
                 threshold=0.5,
                 include_g0 = FALSE,
                 do_sctransform=TRUE,
                 assay='SCT',
                 species='human',
                 gene_id='ensembl',
                 spatial = FALSE)

```
Standard output:  
```r
Running ccAFv2:
  Redoing SCTransform to ensure maximum overlap with classifier genes...
  Total possible marker genes for this classifier: 861
    Marker genes present in this dataset: 854
    Missing marker genes in this dataset: 7
  Predicting cell cycle state probabilities...
217/217 [==============================] - 1s 4ms/step
  Choosing cell cycle state...
  Adding probabilities and predictions to metadata
Done
```

# Cell Cycle Classification Results

The results of the cell cycle classification are stored in the Seurat
object metadata. The likelihoods for each cell cycle state can be found
with the labels of each cell cycle state ( 'Neural.G0', 'G1', 'Late.G1', 'S',
'S.G2', 'G2.M', and 'M.Early.G1') and the classification for each cell
can be found in 'ccAFv2'.

Note that when 'include_g0 = FALSE' there will still be a calculated likelihood for 'Neural.G0', however none of the predictions in the 'ccAFv2' column will be 'G0'. 
Check to ensure 'ccAFv2' colum exists:
```r
colnames(seurat_obj@meta.data)
```
```r
 [1] "orig.ident"   "nCount_RNA"   "nFeature_RNA" "percent.mito" "nCount_SCT"
 [6] "nFeature_SCT" "Neural.G0"    "G1"           "G2.M"         "Late.G1"
[11] "M.Early.G1"   "S"            "S.G2"         "ccAFv2"
```

Show the first few rows of meta.data to see the likelihood estimates for each state, and ccAFv2 predctions in 'ccAFv2' column:
```r
head(seurat_obj@meta.data)
```
```r
                      orig.ident nCount_RNA nFeature_RNA percent.mito
AAACCCAAGTCAGAGC-1 SeuratProject      21363         4216   0.03318822
AAACCCACAGCCTACG-1 SeuratProject      15637         3660   0.04137622
AAACCCACATTACTCT-1 SeuratProject       8409         2637   0.02628137
AAACCCAGTACGATGG-1 SeuratProject       9909         2936   0.03380765
AAACCCAGTGTCTCCT-1 SeuratProject      27986         5085   0.02797827
AAACCCATCCGTGTAA-1 SeuratProject      27078         5055   0.02156732
                   nCount_SCT nFeature_SCT    Neural.G0          G1
AAACCCAAGTCAGAGC-1      20895         4216 1.483502e-02 0.976622045
AAACCCACAGCCTACG-1      19656         3660 9.019502e-03 0.987188399
AAACCCACATTACTCT-1      19332         3013 7.362963e-05 0.001569535
AAACCCAGTACGATGG-1      19859         3153 6.343876e-04 0.009004775
AAACCCAGTGTCTCCT-1      21537         5066 9.892584e-03 0.988429546
AAACCCATCCGTGTAA-1      21558         5044 4.407989e-03 0.982389688
                           G2.M      Late.G1   M.Early.G1            S
AAACCCAAGTCAGAGC-1 6.585211e-03 6.792923e-04 8.019377e-06 0.0010036188
AAACCCACAGCCTACG-1 2.166949e-03 1.017885e-03 4.980507e-10 0.0005518190
AAACCCACATTACTCT-1 5.621322e-05 2.384315e-04 9.968277e-01 0.0011621476
AAACCCAGTACGATGG-1 7.655303e-03 1.499089e-05 9.810986e-01 0.0013027426
AAACCCAGTGTCTCCT-1 9.588006e-04 1.518984e-04 1.379210e-08 0.0005500166
AAACCCATCCGTGTAA-1 4.972871e-03 6.032872e-03 5.351694e-07 0.0021027967
                           S.G2     ccAFv2
AAACCCAAGTCAGAGC-1 2.667967e-04         G1
AAACCCACAGCCTACG-1 5.545938e-05         G1
AAACCCACATTACTCT-1 7.239843e-05 M/Early G1
AAACCCAGTACGATGG-1 2.890851e-04 M/Early G1
AAACCCAGTGTCTCCT-1 1.724525e-05         G1
AAACCCATCCGTGTAA-1 9.308949e-05         G1
```

Count the number of cells in each cell cycle state:
```r
table(seurat_obj@meta.data$ccAFv2)
```
```r
 Neural G0         G1    Late G1          S       S/G2       G2/M M/Early G1
        31       4495        133        365        580        687        916
```


# Setup to plot a UMAP
This data object doesn't yet have UMAP dimensionality reduction applied to it. We will do that using the code below:
```r
seurat_obj = RunPCA(seurat_obj)
seurat_obj = FindNeighbors(seurat_obj)
seurat_obj = RunUMAP(seurat_obj, dims=1:10)
```


# Plotting Cell Cycle States
Plot a UMAP with cell cycle states
```r
DimPlot.ccAFv2(seurat_obj)
```
![Dimplot_1]({{ site.baseurl }}/images/HSkMSCs.png)

