---
layout: default
title: DimPlot.ccAFv2
parent: API
nav_order: 5
---
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
