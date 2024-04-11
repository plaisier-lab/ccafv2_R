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
pdf('ccAFv2_Thresholds_V2.pdf')
ThresholdPlot(seurat_obj)
dev.off()



library(ggplot2)

# 
predictions1 = seurat_obj@meta.data[,c('Neural.G0','G1','Late.G1','S','S.G2','G2.M','M.Early.G1')]
CellCycleState = data.frame(factor(colnames(predictions1)[apply(predictions1,1,which.max)], levels=c('Neural.G0','G1','Late.G1','S','S.G2','G2.M','M.Early.G1','Unknown')), row.names = rownames(predictions1))
colnames(CellCycleState) = 'ccAFv2'
dfall = data.frame(table(CellCycleState)/nrow(CellCycleState))
dfall[,'threshold'] = 0
for(cutoff in c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)) {
    CellCycleState = data.frame(factor(colnames(predictions1)[apply(predictions1,1,which.max)], levels=c('Neural.G0','G1','Late.G1','S','S.G2','G2.M','M.Early.G1','Unknown')), row.names = rownames(predictions1))
    colnames(CellCycleState) = 'ccAFv2'
    CellCycleState[which(apply(predictions1,1,max)<cutoff),'ccAFv2'] = 'Unknown'
    df1 = data.frame(table(CellCycleState)/nrow(CellCycleState))
    df1[,'Threshold'] = as.character(cutoff)
    dfall = rbind(dfall, df1)
}

pdf('ccAFv2_Thresholds.pdf')
p1 = ggplot(dfall) + geom_bar(aes(x = threshold, y = Freq, fill = ccAFv2), position = "stack", stat = "identity") + scale_fill_manual(values = c('G1' = '#f37f73', 'G2.M' = '#3db270', 'Late.G1' = '#1fb1a9', 'M.Early.G1' = '#6d90ca', 'Neural.G0' = '#d9a428', 'S' = '#8571b2', 'S.G2' = '#db7092', 'Unknown' = '#CCCCCC')) + theme_minimal()
print(p1)
dev.off()

