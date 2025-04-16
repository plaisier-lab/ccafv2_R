---
layout: default
title: ThresholdPlot
parent: API
nav_order: 2
---
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
| seurat_obj | NA              | A Seurat object containing metadata columns with cell cycle prediction scores for the following states: 'Neural.G0', 'G1', 'Late.G1', 'S', 'S.G2', 'G2.M', 'M.Early.G1' |
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
Use](https://plaisier-lab.github.io/ccafv2_R/src/Choosing_Threshold.html#selecting-a-likelihood-threshold)
