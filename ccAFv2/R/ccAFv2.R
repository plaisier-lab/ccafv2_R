#' Predict Cell Cycle
#'
#' This function predicts the cell cycle state for each cell in the object
#' using the ccAFv2 cell cycle classifier. [Which cell cycle states/phases]
#' [How to interpret data]
#'
#' @param Seurat object that should have ccAFv2 cell cycle states predicted.
#' @return Seurat object with ccAFv2 calls and probabilities for each cell cycle state.
#' @export
PredictCellCycle = function(seurat1, cutoff=0.5) {
    print('Running ccAFv2:')
    # Load model and marker genes
    ccAFv2 = keras:::load_model_hdf5(system.file('extdata', 'ccAFv2_model.h5', package='ccAFv2'))
    mgenes = read.csv(system.file('extdata', 'ccAFv2_genes.csv', package='ccAFv2'))[,2]
    
    # Subset data marker genes to marker genes included in classification
    seurat_subset = Seurat::subset(seurat1, features = mgenes)

    # Find missing genes and assign 0s to each cell
    print(paste0('  Marker genes for this classifier: ', length(mgenes)))
    missing_genes = setdiff(mgenes, rownames(seurat_subset[['SCT']]@data))
    print(paste0('  Marker genes present in this dataset: ', nrow(seurat_subset[['SCT']]@data)))
    print(paste0('  Missing marker genes in this dataset: ', length(missing_genes)))
    ## Add ERROR later to warn if not enough marker genes ##
    tmp = matrix(0,nrow=length(missing_genes), ncol=ncol(seurat1))
    rownames(tmp) = missing_genes
    colnames(tmp) = colnames(seurat_subset)
    input_mat = seurat_subset[['SCT']]@data
    input_mat_scaled = t(scale(t(as.matrix(input_mat))))
    input_mat_scaled_add_missing_genes = rbind(input_mat_scaled, tmp)[mgenes,]
    print(paste0('  Predicting cell cycle state probabilities...'))
    predictions1 = predict(ccAFv2, t(input_mat_scaled_add_missing_genes))
    colnames(predictions1) = c('G1', 'G1/other', 'G2/M', 'Late G1', 'M/Early G1', 'Neural G0', 'S', 'S/G2')
    rownames(predictions1) = colnames(seurat1)
    df1 = data.frame(predictions1)
    print(paste0('  Choosing cell cycle state...'))
    CellCycleState = data.frame(colnames(predictions1)[apply(predictions1,1,which.max)], row.names = rownames(predictions1))
    colnames(CellCycleState) = 'ccAFv2'
    df1[,'ccAFv2'] = CellCycleState$ccAFv2
    df1[which(apply(predictions1,1,max)<cutoff),'ccAFv2'] = 'Unknown'
    print(paste0('  Adding probabilitities and predictions to metadata'))
    seurat1 = Seurat::AddMetaData(object = seurat1, metadata = df1)
    print('Done')
    return(seurat1)
}

#' DimPlot of ccAFv2 predictions with standard colors
#'
#' This function predicts the cell cycle state for each cell in the object
#' using the ccAFv2 cell cycle classifier. [Which cell cycle states/phases]
#' [How to interpret data]
#'
#' @param Seurat object that should have ccAFv2 cell cycle states predicted.
#' @return Seurat object with ccAFv2 calls and probabilities for each cell cycle state.
#' @export
DimPlot.ccAFv2 = function(seurat1, ...) {
    dp1 = Seurat::DimPlot(seurat1, group.by='ccAFv2', cols = c('G1' = '#f37f73', 'G1/other' = '#9aca3c', 'G2/M' = '#3db270', 'Late G1' = '#1fb1a9','M/Early G1' = '#6d90ca', 'Neural G0' = '#d9a428', 'S' = '#8571b2', 'S/G2' = '#db7092'), ...)
    return(dp1)
}
