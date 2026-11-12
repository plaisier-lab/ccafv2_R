library(Seurat)
test_that("projection works", {
  #U5_data = readRDS('data/U5_normalized_ensembl.rds')
  U5_training_data = read.csv('data/training_data_U5_7_28_2026.csv', header = T, row.names=1)
  U5_training_data = as(as.matrix(U5_training_data), "dgCMatrix")
  U5_data = CreateSeuratObject(counts = U5_training_data, project = 'U5', assay = "SCT")
  #U5_data = SetAssayData(U5_data, slot = "scale.data", new.data = U5_training_data)
  #U5_data[['SCT']]@scale.data = as.matrix(U5_training_data)
  LayerData(U5_data, assay='SCT', layer='scale.data') = U5_training_data
  U5_ref = readRDS(system.file("extdata", "seurat_ref_U5_hNSC_SCT_7_29_2026.rds", package = "ccAFv2"))
  #U5_data = SCTransform(U5_data, return.only.var.genes=FALSE, verbose=FALSE)
  U5_data2 = suppressMessages(suppressWarnings(ProjectCycleFromSeurat(
    U5_data,
    U5_ref,
    assay = 'SCT',
    layer = 'scale.data',
    gene_col = 1,
    reduction_name = 'cc_Projection',
    reduction_key = 'CYCLE_',
    position_col = 'thetaPos',
    center.pc1 = 0,
    center.pc2 = 0
  )))
  expect_equal(U5_data2@meta.data$thetaPos[length(U5_data2@meta.data$thetaPos)], 3.1139751306)
})
