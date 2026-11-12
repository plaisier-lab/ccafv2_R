test_that("scaled values of bins using vonMises is correct", {
  U5_data = readRDS('data/U5_normalized_ensembl.rds')
  U5_data = suppressMessages(suppressWarnings(SCTransform(U5_data, return.only.var.genes = F, verbose = F)))
  U5_ref = readRDS(system.file("extdata", "seurat_ref_U5_hNSC_SCT_7_29_2026.rds", package = "ccAFv2"))
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
  binProbs = AnglesToVonMisesBins(U5_data2@meta.data$thetaPos)
  scaledBinProbs = t(.scale(t(binProbs)))
  #print(head(binProbs))
  expect_equal(as.numeric(scaledBinProbs[1,]), c(1.091647426351784, 0.15612371667585365, -0.5052900829024036, -0.5216521728838878, -0.541440590956184, -0.5957156022368719, -0.5103879095091882, -0.39465152018566396, 0.33372008955700233, 1.1564861921070022))
})
