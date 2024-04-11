# docker run -it -v '/media/old_home/home/cplaisier/ccAFv2_test:/files' cplaisier/ccafv2_extra

devtools::install_github("plaisier-lab/ccafv2_R/ccAFv2")

library(ccAFv2)
library(Seurat)

reticulate:::use_python('/usr/bin/python3')

# Load Data
setwd('/files')
seurat_obj = readRDS('U5_normalized_ensembl.rds')

# Predict cell cycle states
seurat_obj = PredictCellCycle(seurat_obj)

# Plot DimPlot colorized by cell cycle states
pdf('ccAFv2_DimPlot.pdf')
DimPlot.ccAFv2(seurat_obj)
dev.off()

# Plot DimPlot colorized by cell cycle states
pdf('ccAFv2_ThresholdPlot.pdf')
ThresholdPlot(seurat_obj)
dev.off()

# Regress out cell cycle
seurat_obj = PrepareForCellCycleRegression(seurat_obj)
seurat_obj = SCTransform(data2, vars.to.regress = c("Late.G1_exprs1", "S_exprs2", "S.G2_exprs3", "G2.M_exprs4", "M.Early.G1_exprs5"))

# Plot cell cycle regressed UMAP
seurat_obj = RunPCA(seurat_obj)
seurat_obj = RunUMAP(suerat_obj)
pdf('ccAFv2_DimPlot_regressed.pdf')
DimPlot.ccAFv2(seurat_obj)
dev.off()


