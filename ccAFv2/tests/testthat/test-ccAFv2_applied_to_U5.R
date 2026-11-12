library(Seurat)
library(ggplot2)
library(pROC)
test_that("running ccAFv2 function", {
  seurat_obj = readRDS('data/U5_normalized_ensembl.rds')
  test_cells = as.character(read.csv('data/test_hold_out_cells_19.csv', header=T)[,2])
  out1 = capture.output({seurat_obj = suppressMessages(suppressWarnings(PredictCellCycle(seurat_obj, include_g0 = T)))}, type='output')
  convertMe = c('Neural.G0', 'G1', 'Late.G1', 'S', 'S.G2', 'G2.M', 'M.Early.G1')
  names(convertMe) = c('Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1')
  tmp_ccAF = seurat_obj@meta.data[test_cells,'ccAF']
  tmp_ccAF = convertMe[tmp_ccAF]
  tmp_ccAF = factor(tmp_ccAF, levels = c('Neural.G0', 'G1', 'Late.G1', 'S', 'S.G2', 'G2.M', 'M.Early.G1'))
  #print(table(seurat_obj@meta.data[test_cells,'ccAFv2'], tmp_ccAF))
  auc_per_class = sapply(c('Neural.G0', 'G1', 'Late.G1', 'S', 'S.G2', 'G2.M', 'M.Early.G1'), function(class_label) {
    # 1. Create binary response: 1 if the true class is the current class, 0 otherwise
    binary_true = ifelse(tmp_ccAF == class_label, 1, 0)
    
    # 2. Extract probabilities for the current class
    binary_probs = seurat_obj@meta.data[test_cells, class_label]
    
    # 3. Compute ROC and AUC for this specific class
    roc_obj = roc(binary_true, binary_probs, quiet = TRUE)
    return(auc(roc_obj))
  })
  #print(auc_per_class)
  #roc_obj = multiclass.roc(predictor = seurat_obj@meta.data[test_cells,c('Neural.G0', 'G1', 'Late.G1', 'S', 'S.G2', 'G2.M', 'M.Early.G1')], response = tmp_ccAF)
  #print(auc(roc_obj))
  tmp1 = as.character(head(seurat_obj@meta.data[test_cells,'ccAF'], n=10))
  tmp1[2] = 'G1'
  expect_gte(min(auc_per_class), 0.9)
  expect_equal(as.character(head(seurat_obj@meta.data[test_cells,'ccAFv2'], n=10)), tmp1)
})
