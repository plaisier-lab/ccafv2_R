---
layout: default
title: PredictCellCycle
parent: API
nav_order: 1
---

## PredictCellCycle

This function predicts the cell cycle state for each cell in the object
using the ccAFv2 cell cycle classifier. The possible cell cycle states
that ccAFv2 can predict are: Neural G0, G1, G1/other, Late G1, S, S/G2,
G2/M, and M/Early G1.

``` r        
PredictCellCycle(
  seurat_obj,
  threshold = 0.5,
  include_g0 = FALSE,
  do_sctransform = TRUE,
  assay = "SCT",
  species = "human",
  gene_id = "ensembl",
  spatial = FALSE
 )
```

### Arguments

| Arguments      | Range of Values        | Description                                                                                           |
|-------------------|-------------------|-----------------------------------|
| seurat_obj     |           NA             | a seurat object must be supplied to classify, no default                                              |
| threshold      | Numeric (0 to 1)       | the value used to threshold the likelihoods, default is 0.5                                           |
| include_g0     | TRUE or FALSE          | whether to provide Neural G0 calls, or collapse G1 and Neural G0 into G0/G1, default is FALSE         |
| do_sctransform | TRUE or FALSE          | whether to do SCTransform before classifying, default is TRUE                                         |
| assay          | 'SCT', 'RNA', etc.     | which seurat_obj assay to use for classification, helpful if data is pre-normalized, default is 'SCT' |
| species        | 'human' or 'mouse',    | from which species did the samples originate, default is 'human'                                      |
| gene_id        | 'ensembl' or 'symbol', | what type of gene ID is used, default is 'ensembl'                                                    |
| spatial        | TRUE or FALSE          | whether the data is spatial, default is FALSE                                                         |

### Value

Returns the input Seurat object, with:

-   Cell cycle prediction probabilities added to the metadata (one
    column per class).

-   A 'ccAFv2' column in metadata containing the assigned cell cycle
    state per cell.

### Example

``` r       
# Run classifier and add predictions to Seurat metadata
seurat_obj = PredictCellCycle(seurat_obj)
```

[Tutorial
example](https://plaisier-lab.github.io/ccafv2_R/src/U5.html#marker-genes)
