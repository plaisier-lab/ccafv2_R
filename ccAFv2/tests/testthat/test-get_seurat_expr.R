test_that("grabbing seurat expression works", {
  U5_data = readRDS('data/U5_normalized_ensembl.rds')
  U5_data = suppressMessages(suppressWarnings(SCTransform(U5_data, return.only.var.genes = F, verbose = F)))
  exprs = GetSeuratExpr(U5_data, assay='SCT')
  expect_equal(dim(exprs), c(13866, 2962))
})
