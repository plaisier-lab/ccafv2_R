---
layout: default
title: API
nav_order: 4
has_children: true
---

# API
### Table of Contents

- **Prediction & Thresholding**
  - [`PredictCellCycle`](#predictcellcycle)  
    Run the ccAFv2 classifier and assign cell cycle states.
  - [`AdjustCellCycleThreshold`](#adjustcellcyclethreshold)  
    Re-assign states based on a new prediction confidence threshold.
  - [`ThresholdPlot`](#thresholdplot)  
    Visualize how prediction thresholds affect state assignment.

- **Visualization**
  - [`DimPlot.ccAFv2`](#dimplotccafv2)  
    Plot predicted cell cycle states on a UMAP or PCA.
  - [`SpatialDimPlot.ccAFv2`](#spatialdimplotccafv2)  
    Visualize predicted states in spatial transcriptomics data.

- **Data Preparation**
  - [`PrepareForCellCycleRegression`](#prepareforcellcycleregression)  
    Compute module scores from marker gene sets for regression.
