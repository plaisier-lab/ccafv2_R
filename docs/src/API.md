---
layout: default
title: API
nav_order: 4
has_children: true
toc: false
---

# API
### Table of Contents

- **Prediction & Thresholding**
  - [`PredictCellCycle`](https://plaisier-lab.github.io/ccafv2_R/src/PredictCellCycle.html) 
    Run the ccAFv2 classifier and assign cell cycle states.
  - [`AdjustCellCycleThreshold`](https://plaisier-lab.github.io/ccafv2_R/src/adjustcellcyclethreshold)  
    Re-assign states based on a new prediction confidence threshold.
  - [`ThresholdPlot`](https://plaisier-lab.github.io/ccafv2_R/src/thresholdplot)  
    Visualize how prediction thresholds affect state assignment.

- **Visualization**
  - [`DimPlot.ccAFv2`](https://plaisier-lab.github.io/ccafv2_R/src/dimplotccafv2)  
    Plot predicted cell cycle states on a UMAP or PCA.
  - [`SpatialDimPlot.ccAFv2`](https://plaisier-lab.github.io/ccafv2_R/src/spatialdimplotccafv2)  
    Visualize predicted states in spatial transcriptomics data.

- **Data Preparation**
  - [`PrepareForCellCycleRegression`](https://plaisier-lab.github.io/ccafv2_R/src/prepareforcellcycleregression)  
    Compute module scores from marker gene sets for regression.
