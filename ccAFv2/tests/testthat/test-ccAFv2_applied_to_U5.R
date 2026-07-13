library(Seurat)
library(ggplot2)
test_that("scaled values of bins using vonMises is correct", {
  seurat_obj = readRDS('data/U5_normalized_ensembl.rds')
  test_cells = as.character(read.csv('data/test_hold_out_cells_8.csv', header=T)[,'CellID'])
  seurat_obj = suppressMessages(suppressWarnings(PredictCellCycle(seurat_obj[,test_cells], include_g0=TRUE)))
  #print(table(seurat_obj@meta.data$ccAFv2, seurat_obj@meta.data$ccAF))
  expect_equal(as.character(head(seurat_obj@meta.data$ccAFv2, n=10)), c('Late G1', 'G1', 'S/G2', 'M/Early G1', 'S', 'M/Early G1', 'G1', 'G2/M', 'M/Early G1', 'G1'))
})
