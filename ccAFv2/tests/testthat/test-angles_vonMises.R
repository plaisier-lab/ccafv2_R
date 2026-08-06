test_that("binning using vonMises is correct", {
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
  #print(head(binProbs))

  expect_equal(as.numeric(binProbs[1,]), c(0.4527713397868222, 0.22070704747892222, 0.03406083720417943, 0.0029953642492716055, 0.00037839667396732967, 0.00017355269458152478, 0.000428537146549113 ,0.003617001300473712, 0.04058066689418351, 0.24428725657104935))
})
