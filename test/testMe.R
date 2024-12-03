# docker run -it -v '/media/old_home/home/cplaisier/ccAFv2_test:/files' cplaisier/ccafv2_seurat4

devtools::install_github("plaisier-lab/ccafv2_R/ccAFv2")

library(ccAFv2)
library(Seurat)

reticulate:::use_python('/usr/bin/python3')

### Single cell or nuclei

# Load Data
setwd('/files')
seurat_obj = readRDS('U5_normalized_ensembl.rds')


##################
## No Neural G0 ##
##################

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

# Adjust ccAFv2 threshold to 0.9
seurat_obj = AdjustCellCycleThreshold(seurat_obj, threshold=0.9)

# Plot DimPlot colorized by cell cycle states
pdf('ccAFv2_DimPlot_T_0.9.pdf')
DimPlot.ccAFv2(seurat_obj)
dev.off()

# Adjust ccAFv2 threshold back to 0.5
seurat_obj = AdjustCellCycleThreshold(seurat_obj, threshold=0.5)

# Plot DimPlot colorized by cell cycle states
pdf('ccAFv2_DimPlot_T_0.5.pdf')
DimPlot.ccAFv2(seurat_obj)
dev.off()


####################
## With Neural G0 ##
####################

# Predict cell cycle states
seurat_obj = PredictCellCycle(seurat_obj, include_g0=TRUE)

# Plot DimPlot colorized by cell cycle states
pdf('ccAFv2_DimPlot_wNeuralG0.pdf')
DimPlot.ccAFv2(seurat_obj)
dev.off()

# Plot DimPlot colorized by cell cycle states
pdf('ccAFv2_ThresholdPlot_wNeuralG0.pdf')
ThresholdPlot(seurat_obj)
dev.off()

# Adjust ccAFv2 threshold to 0.9
seurat_obj = AdjustCellCycleThreshold(seurat_obj, threshold=0.9, include_g0=TRUE)

# Plot DimPlot colorized by cell cycle states
pdf('ccAFv2_DimPlot_T_0.9_wNeuralG0.pdf')
DimPlot.ccAFv2(seurat_obj)
dev.off()

# Adjust ccAFv2 threshold back to 0.5
seurat_obj = AdjustCellCycleThreshold(seurat_obj, threshold=0.5, include_g0=TRUE)

# Plot DimPlot colorized by cell cycle states
pdf('ccAFv2_DimPlot_T_0.5_wNeuralG0.pdf')
DimPlot.ccAFv2(seurat_obj)
dev.off()


###########################
## Cell cycle regression ##
###########################


# Regress out cell cycle
seurat_obj = PrepareForCellCycleRegression(seurat_obj)
seurat_obj = SCTransform(seurat_obj, vars.to.regress = c("Late.G1_exprs1", "S_exprs2", "S.G2_exprs3", "G2.M_exprs4", "M.Early.G1_exprs5"))

# Plot cell cycle regressed UMAP
seurat_obj = RunPCA(seurat_obj)
seurat_obj = RunUMAP(seurat_obj, dims=1:10)
pdf('ccAFv2_DimPlot_regressed.pdf')
DimPlot.ccAFv2(seurat_obj)
dev.off()


#############
## Spatial ##
#############

# Load Data
setwd('/files')
spatial_obj = readRDS('GSM6736780_Spatial_10x_PCW4_20220122_slice1.rds')

# Predict cell cycle states
spatial_obj = PredictCellCycle(spatial_obj, species='human', gene_id='symbol', spatial=TRUE)

# Plot cell cycle states onto spatial image
pdf('ccAFv2_SpatialDimPlot_slice1.pdf')
SpatialDimPlot.ccAFv2(spatial_obj) + theme(legend.position = "right")
dev.off()


