---
layout: default
title: Cell Cycle Regression 
nav_order: 10

---
# Characterizing Grade 3 Oligodendroglioma
We characterized the transcriptomes of cells from a grade 3 oligodendroglioma (Garcia et al. 2026) cultured in vitro using our standrd single cell pipeline: [predict the cell cycle]([https://plaisier-lab.github.io/ccafv2_R/src/single.html)

```r
# Load libraries
library(ccAFv2)
library(Seurat)

# Set working directory
setwd('/files/U5_test')

# Download data file
download.file(url="https://zenodo.org/records/24207612/files/BT237_GF_normalized_ensembl.rds?download=1", destfile="BT237_GF_normalized_ensembl.rds", mode="wb")

# Load up sample data
seurat_obj = readRDS('BT237_GF_normalized_ensembl.rds')

# Run ccAFv2
seurat_obj = PredictCellCycle(seurat_obj,
                 threshold=0,
                 include_g0 = TRUE,
                 do_sctransform=TRUE,
                 assay='SCT',
                 layer = 'scale.data',
                 species='human',
                 gene_id='ensembl',
                 spatial = FALSE)

# Check to ensure 'ccAFv2' column exists
colnames(seurat_obj@meta.data)

# Show the first few rows of meta.data to see the likelihood estimates for each
# state, and ccAFv2 predictions in the 'ccAFv2' column
head(seurat_obj@meta.data)

# Number of cells in each cell cycle state
table(seurat_obj@meta.data$ccAFv2)

# Setup to plot a UMAP
seurat_obj = RunPCA(seurat_obj)
seurat_obj = FindNeighbors(seurat_obj)
seurat_obj = RunUMAP(seurat_obj, dims=1:10)

# Plot the UMAP with ccAFv2 predicitions
png('BT237_GF_UMAP.png', width = 400, height = 400)
DimPlot.ccAFv2(seurat_obj)
dev.off()
```
The resulting UMAP shows a significant cell cycle bias:
![PreRegression_Dimplot]({{ site.baseurl }}/images/BT237_GF_UMAP.png)

# Cell Cycle Regression
Since the cell cycle significantly impacts gene expression, it is common practice to remove its effects and use the residual variance for further analysis. We support this using ccAFv2 marker genes, starting by calculating expression module scores for G0, Late G1, S, S/G2, G2/M, and M/Early G1 phases.

## Prepare For Cell Cycle Regression
```r
seurat_obj = PrepareForCellCycleRegression(seurat_obj)
```

## Conduct Cell Cycle Regression
```
seurat_obj = SCTransform(seurat_obj, vars.to.regress = c("Neural.G0","Late.G1_exprs1", "S_exprs2", "S.G2_exprs3", "G2.M_exprs4", "M.Early.G1_exprs5"))
```

## See the Impact of Cell Cycle Regression
Removing the cell cycle from the BT237 tumor cells leads to an intersting multiprong distribution.

### Redo the UMAP After Regression
It is important to the redo the UMAP after regressing out the cell cycle to see what impact the cell cycle regression has had on the data.
```r
seurat_obj = RunPCA(seurat_obj)
seurat_obj = RunUMAP(seurat_obj, dims=1:10)
```

### Plot After Regression
```r
DimPlot.ccAFv2(seurat_obj)
```

![AfterRegression]({{ site.baseurl }}/images/BT237_umap_regressOut.png)
