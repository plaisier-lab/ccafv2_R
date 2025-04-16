---
layout: default
title: AdjustCellCycleThreshold
parent: API
nav_order: 3
---

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

Returns the modified Seurat object, updating the 'ccAFv2' metadata
column to reflect cell cycle state predictions at the new threshold.

### Example

``` r
# Adjust prediction stringency
AdjustCellCycleThreshold(seurat_obj)
```

[Example
Use](https://plaisier-lab.github.io/ccafv2_R/src/Choosing_Threshold.html#selecting-a-likelihood-threshold)
