---
layout: default
title: PrepareForCellCycleRegression
parent: API
nav_order: 4
---

## PrepareForCellCycleRegression

This function calculates module scores for a curated set of cell
cycle-related gene clusters from the `ccAFv2` gene set. These scores can
be used to regress out cell cycle effects from single-cell expression
data.

``` r        
PrepareForCellCycleRegression(
  seurat_obj,
  assay = "SCT",
  species = "human",
  gene_id = "ensembl"
)
```

| Arguments  | Range of Values       | Description                                                                                           |
|-------------------|------------------|------------------------------------|
| seurat_obj | NA                    | a seurat object must be supplied to classify, no default                                              |
| assay      | 'SCT','RNA', etc.     | which seurat_obj assay to use for classification, helpful if data is pre-normalized, default is 'SCT' |
| species    | 'human' or 'mouse'    | from which species did the samples originate, default is 'human'                                      |
| gene_id    | 'ensembl' or 'symbol' | what type of gene ID is used, default is 'ensembl'                                                    |

### Value

Returns the input Seurat object, with five new metadata columns
containing module scores for the following cell cycle clusters (Late.G1,
S, S.G2, G2.M, M.Early.G1)

These scores represent average expression of marker genes for each phase
and can be used for downstream cell cycle regression.

### Example

``` r
# Prepare Seurat object for cell cycle regression 
PrepareForCellCycleRegression(seurat_obj)
```

[Example
Use](https://plaisier-lab.github.io/ccafv2_R/src/regress.html#cell-cycle-regression)
