test_that("projection works", {
  U5_data = readRDS('data/U5_normalized_ensembl.rds')
  U5_data2 = suppressMessages(suppressWarnings(ProjectCycleFromSeurat(
    U5_data,
    system.file("extdata", "seurat_ref_U5_hNSC.rds", package = "ccAFv2"),
    assay = 'RNA',
    layer = 'scale.data',
    gene_col = 1,
    reduction_name = 'cc_Projection',
    reduction_key = 'CYCLE_',
    position_col = 'thetaPos',
    center.pc1 = 0,
    center.pc2 = 0
  )))
  expect_equal(U5_data2@meta.data$thetaPos[length(U5_data2@meta.data$thetaPos)], 2.6121559229)
})
