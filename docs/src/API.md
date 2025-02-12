---
layout: default
title: API
nav_order: 4
---

# API

## PredictCellCycle

This function predicts the cell cycle state for each cell in the object
using the ccAFv2 cell cycle classifier. The possible cell cycle states
that ccAFv2 can predict are: Neural G0, G1, G1/other, Late G1, S, S/G2,
G2/M, and M/Early G1.

```         
PredictCellCycle(
  seurat_obj,
  threshold = 0.5,
  include_g0 = FALSE,
  do_sctransform = TRUE,
  assay = "SCT",
  species = "human",
  gene_id = "ensembl",
  spatial = FALSE)
```


| Arguments      | Description                                                                                                                 |
|------------------------------------|------------------------------------|
| seurat_obj     | a seurat object must be supplied to classify, no default                                                                    |
| threshold      | the value used to threshold the likelihoods, default is 0.5                                                                 |
| include_g0     | whether to provide Neural G0 calls, or collapse G1 and Neural G0 into G0/G1 (FALSE collapses, TRUE provides Neural G0 calls |
| do_sctransform | whether to do SCTransform before classifying, default is TRUE                                                               |
| assay          | which seurat_obj assay to use for classification, helpful if data is pre-normalized, default is 'SCT'                       |
| species        | from which species did the samples originate, either 'human' or 'mouse', default is 'human'                                 |
| gene_id        | what type of gene ID is used, either 'ensembl' or 'symbol', defaults to 'ensembl'                                           |
| spatial        | whether the data is spatial, defaults to FALSE                                                                              |

[Example Use](https://rlhoove1.github.io/tryingmybest/src/single.html#marker-genes)

## ThresholdPlot

This function plots the distribution of cell cycle predictions for a
range of thresholds as a barplot colorized using the standard cell cycle
state colors.

```         
ThresholdPlot(seurat_obj, ...)
```


| Argument   | Description                                               |
|------------|-----------------------------------------------------------|
| seurat_obj | A Seurat object must be supplied to classify, no default. |

[Example
Use](https://rlhoove1.github.io/tryingmybest/src/Choosing_Threshold.html)

## AdjustCellCycleThreshold

This function allows users to adjust the threshold applied to ccAFv2
predictions. The user can utilize the ThresholdPlot to see the effect of
increasing threshold values have on the number of 'Unknown' cell calls.

```         
AdjustCellCycleThreshold(seurat_obj, threshold = 0.5, include_g0 = FALSE)
```

| Arguments  | Description                                                                                                                 |
|------------------------------------|------------------------------------|
| seurat_obj | a seurat object must be supplied to classify, no default                                                                    |
| threshold  | the value used to threshold the likelihoods, default is 0.5                                                                 |
| include_g0 | whether to provide Neural G0 calls, or collapse G1 and Neural G0 into G0/G1 (FALSE collapses, TRUE provides Neural G0 calls |

[Example
Use](https://rlhoove1.github.io/tryingmybest/src/Choosing_Threshold.html)

## PrepareForCellCycleRegression

This function computes module scores for each cell cycle state

```         
PrepareForCellCycleRegression(
  seurat_obj,
  assay = "SCT",
  species = "human",
  gene_id = "ensembl"
)
```

| Arguments  | Description                                                                                           |
|------------------------------------|------------------------------------|
| seurat_obj | a seurat object must be supplied to classify, no default                                              |
| assay      | which seurat_obj assay to use for classification, helpful if data is pre-normalized, default is 'SCT' |
| species    | from which species did the samples originate, either 'human' or 'mouse', default is 'human'           |
| gene_id    | what type of gene ID is used, either 'ensembl' or 'symbol', defaults to 'ensembl'                     |

[Example Use](https://rlhoove1.github.io/tryingmybest/src/regress.html#cell-cycle-regression)

## DimPlot.ccAFv2

This function plots the cell cycle onto a DimPlot for single cell or
nuclei data.

```         
DimPlot.ccAFv2(seurat_obj, ...)
```

| Arguments  | Description                                                                                          |
|------------------------------------|------------------------------------|
| seurat_obj | a seurat object must be supplied to classify, no default                                             |
| ...        | This function supports all [**DimPlot**](https://satijalab.org/seurat/reference/dimplot) parameters. |

[Example Use](https://rlhoove1.github.io/tryingmybest/src/single.html#plotting-cell-cycle-states)

## SpatialDimPlot.ccAFv2

This function plots the cell cycle onto a DimPlot for spatial data.

```         
SpatialDimPlot.ccAFv2(seurat_obj, ...)
```

| Arguments  | Description                                                                                                                                      |
|------------------------------------|------------------------------------|
| seurat_obj | a seurat object must be supplied to classify, no default                                                                                         |
| ...        | This function supports all [**SpatialDimPlot**](https://satijalab.org/seurat/reference/spatialplot) parameters except **group.by** and **cols**. |

[Example Use](https://rlhoove1.github.io/tryingmybest/src/spatial.html#plotting-cell-cycle-states-onto-images)
