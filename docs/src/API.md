---
layout: default
title: API
nav_order: 4
has_children: true
---

# API

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

-   A `ccAFv2` column in metadata containing the assigned cell cycle
    state per cell.

### Example

``` r       
# Run classifier and add predictions to Seurat metadata
seurat_obj = PredictCellCycle(seurat_obj)
```

[Tutorial
example](https://rlhoove1.github.io/tryingmybest/src/single.html#marker-genes)

## ThresholdPlot

This function visualizes how classification confidence impacts the
assignment of cells to specific cell cycle states. It does so by
plotting the distribution of predicted cell cycle states across a range
of score thresholds, using a stacked barplot with standardized cell
cycle state colors.

At each threshold:

-   Each cell is assigned to the state with the highest prediction
    score, only if that score exceeds the threshold.

-   Cells with no scores above the threshold are labeled as "Unknown".

-   The proportion of cells in each state is then calculated and
    plotted.

``` r        
ThresholdPlot(seurat_obj, ...)
```

| Argument   | Range of Values | Description                                                                                                                                                                           |
|-----------------|---------------|-------------------------------------------------------|
| seurat_obj | NA              | A Seurat object containing metadata columns with cell cycle prediction scores for the following states: `'Neural.G0'`, `'G1'`, `'Late.G1'`, `'S'`, `'S.G2'`, `'G2.M'`, `'M.Early.G1'` |
| ...        | ...             | Additional arguments passed using `ggplot` functions for further customization.                                                                                                       |

### Value

Returns a ggplot object showing the proportion of cells classified into
each cell cycle state at different classification score thresholds.
Cells with a maximum score below a given threshold are labeled as
"Unknown" for that threshold.

### Example

``` r
ThresholdPlot(seurat_obj, ...)
```

[Example
Use](https://rlhoove1.github.io/tryingmybest/src/Choosing_Threshold.html)

## AdjustCellCycleThreshold

This function allows users to adjust the threshold applied to ccAFv2
predictions. The user can utilize the ThresholdPlot to see the effect of
increasing threshold values have on the number of 'Unknown' cell calls.

``` r        
AdjustCellCycleThreshold(seurat_obj, threshold = 0.5, include_g0 = FALSE)
```

| Arguments  | Range of Values | Description                                                                                   |
|--------------|--------------|------------------------------------|
| seurat_obj | NA              | a seurat object must be supplied to classify, no default                                      |
| threshold  | Numeric(0-1)    | the value used to threshold the likelihoods, default is 0.5                                   |
| include_g0 | TRUE or FALSE   | whether to provide Neural G0 calls, or collapse G1 and Neural G0 into G0/G1, default is FALSE |

### Value

Returns the **modified `Seurat` object**, updating the `ccAFv2` metadata
column to reflect cell cycle state predictions at the new threshold.

### Example

``` r
# Adjust prediction stringency
AdjustCellCycleThreshold(seurat_obj)
```

[Example
Use](https://rlhoove1.github.io/tryingmybest/src/Choosing_Threshold.html)

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
Use](https://rlhoove1.github.io/tryingmybest/src/regress.html#cell-cycle-regression)

## DimPlot.ccAFv2

This function creates a dimensional reduction plot (e.g., UMAP or t-SNE)
to visualize cell cycle states predicted in the `ccAFv2` metadata field.
It uses Seurat’s `DimPlot` function, with standardized coloring for each
cell cycle phase.

``` r         
DimPlot.ccAFv2(seurat_obj, ...)
```

| Arguments  | Range of Values | Description                                                                                          |
|-------------------|-------------------|------------------------------------|
| seurat_obj | NA              | a seurat object must be supplied to classify, no default                                             |
| ...        | ...             | This function supports all [**DimPlot**](https://satijalab.org/seurat/reference/dimplot) parameters. |

### Value

Returns a **`ggplot`** object produced by Seurat’s `DimPlot` function.

### Example

``` r
DimPlot.ccAFv2(seurat_obj)
```

[Example
Use](https://rlhoove1.github.io/tryingmybest/src/single.html#plotting-cell-cycle-states)

## SpatialDimPlot.ccAFv2

This function generates a spatial visualization of cell cycle states
using Seurat's `SpatialDimPlot`. It maps the `ccAFv2` annotations onto
the spatial coordinates of the tissue section, using standardized colors
for each predicted cell cycle state.

``` r         
SpatialDimPlot.ccAFv2(seurat_obj, ...)
```

| Arguments  | Range of Values | Description                                                                                                                                                                           |
|---------------|--------------------|--------------------------------------|
| seurat_obj | NA              | a seurat object must be supplied to classify, no default                                                                                                                              |
| ...        | ...             | Additional arguments passed to [**SpatialDimPlot**](https://satijalab.org/seurat/reference/spatialplot), except for `group.by` and `cols`, which are set internally by this function. |

### Value

Returns a **`ggplot` object** produced by Seurat’s `SpatialDimPlot`
function.

### Example

``` r
# Assuming seurat_obj contains the `ccAFv2` metadata
SpatialDimPlot.ccAFv2(seurat_obj)
```

[Example
Use](https://rlhoove1.github.io/tryingmybest/src/spatial.html#plotting-cell-cycle-states-onto-images)
