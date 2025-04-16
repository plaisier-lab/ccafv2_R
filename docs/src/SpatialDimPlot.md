---
layout: default
title: SpatialDimPlot.ccAFv2
parent: API
nav_order: 6
---
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
