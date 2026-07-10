test_that("grabbing seurat expression works", {
  U5_data = readRDS('data/U5_normalized_ensembl.rds')
  exprs = GetSeuratExpr(U5_data, assay='SCT')
  expect_equal(dim(exprs), c(13866, 2962))
})
