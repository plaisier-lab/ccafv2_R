test_that("binning using vonMises is correct", {
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
  #print(head(binProbs))
  expect_equal(as.numeric(binProbs[1,]), c(0.175880859,	0.023646516,	0.002059723,	0.000301588,	0.000179271,	0.000564031,	0.00532218,	0.057128087,	0.293331236,	0.441586508))
})
