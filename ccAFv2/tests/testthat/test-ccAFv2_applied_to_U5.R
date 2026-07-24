library(Seurat)
library(ggplot2)
library(pROC)
test_that("scaled values of bins using vonMises is correct", {
  seurat_obj = readRDS('data/U5_normalized_ensembl.rds')
  test_cells = as.character(read.csv('data/test_hold_out_cells_8.csv', header=T)[,'CellID'])
  #seurat_obj = suppressMessages(suppressWarnings(PredictCellCycle(seurat_obj[,test_cells], include_g0=TRUE)))
  seurat_obj = suppressMessages(suppressWarnings(PredictCellCycle(seurat_obj, include_g0=TRUE)))
  write.csv(seurat_obj@meta.data, 'meta_data.csv')
  print(colnames(seurat_obj@meta.data))
  convertMe = c('Neural.G0', 'G1', 'Late.G1', 'S', 'S.G2', 'G2.M', 'M.Early.G1')
  names(convertMe) = c('Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1')
  tmp_ccAF = seurat_obj@meta.data$ccAF
  tmp_ccAF = convertMe[tmp_ccAF]
  tmp_ccAF = factor(tmp_ccAF, levels = c('Neural.G0', 'G1', 'Late.G1', 'S', 'S.G2', 'G2.M', 'M.Early.G1'))
  print(table(seurat_obj@meta.data$ccAFv2, tmp_ccAF))
  cal_lab8 = read.csv('data/calibration_labels_8.csv', header=F)[,1]
  print(length(cal_lab8))
  print(length(seurat_obj@meta.data$ccAFv2))
  print(table(seurat_obj@meta.data$ccAFv2, cal_lab8[1:297]))
  auc_per_class = sapply(c('Neural.G0', 'G1', 'Late.G1', 'S', 'S.G2', 'G2.M', 'M.Early.G1'), function(class_label) {
    # 1. Create binary response: 1 if the true class is the current class, 0 otherwise
    binary_true = ifelse(tmp_ccAF == class_label, 1, 0)
    
    # 2. Extract probabilities for the current class
    binary_probs = seurat_obj@meta.data[, class_label]
    
    # 3. Compute ROC and AUC for this specific class
    roc_obj = roc(binary_true, binary_probs, quiet = TRUE)
    return(auc(roc_obj))
  })
  print(auc_per_class)
  roc_obj = multiclass.roc(predictor = seurat_obj@meta.data[,c('Neural.G0', 'G1', 'Late.G1', 'S', 'S.G2', 'G2.M', 'M.Early.G1')], response = tmp_ccAF)
  print(auc(roc_obj))
  expect_equal(as.character(head(seurat_obj@meta.data$ccAFv2, n=10)), c('Late G1', 'G1', 'S/G2', 'M/Early G1', 'S', 'M/Early G1', 'G1', 'G2/M', 'M/Early G1', 'G1'))
})
