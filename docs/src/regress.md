---
layout: default
title: Cell Cycle Regression 
nav_order: 10

---
# Cell Cycle Regression
The cell cycle strongly influences gene expression, so it's common to
regress out cell cycle effects and use residual variance for further
analysis. We support this using ccAFv2 marker genes, starting by
calculating expression module scores for Late G1, S, S/G2, G2/M, and
M/Early G1 phases.

```r
#Collect expression module scores for the cell cycle states 
seurat_obj = PrepareForCellCycleRegression(seurat_obj)

#Regress these signatures out of the expression data
seurat_obj = SCTransform(seurat_obj, vars.to.regress = c("Late.G1_exprs1", "S_exprs2", "S.G2_exprs3", "G2.M_exprs4", "M.Early.G1_exprs5"))
```
Removing the cell cycle from the U5 hNSCs leads to a random distribution
because the cell cycle is the primary biological signal in this in
vitro-grown cell line.
```r
seurat_obj = RunPCA(seurat_obj)
seurat_obj = RunUMAP(seurat_obj, dims=1:10)
DimPlot.ccAFv2(seurat_obj)
```
![Dimplot_2]({{ site.baseurl }}/images/DimPlot_2.jpeg)
