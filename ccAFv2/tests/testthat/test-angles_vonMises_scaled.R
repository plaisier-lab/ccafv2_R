test_that("scaled values of bins using vonMises is correct", {
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
  binProbs = AnglesToVonMisesBins(U5_data2@meta.data$thetaPos)
  scaledBinProbs = t(.scale(t(binProbs)))
  #print(head(binProbs))
  expect_equal(as.numeric(scaledBinProbs[1,]), c(-0.04708347, -0.588508922, -0.552518564, -0.586031215, -0.640347976, -0.5703681, -0.412050272, 0.568508933, 1.669709428, 1.165504172))
})
