#' Predict Cell Cycle
#'
#' This function predicts the cell cycle state for each cell in the object
#' using the ccAFv2 cell cycle classifier. The possible cell cycle states that
#' ccAFv2 can predict are: Neural G0, G1, G1/other, Late G1, S, S/G2, G2/M,
#' and M/Early G1.
#'
#' The ccAFv2 predicts the cell cycle state for each cell in the object by
#' selecting the cell cycle state for each cell with the maximum cell cycle
#' state probability. If the cell cycle state probability for a cell does not
#' meet the probability cutoff, the cell will receive an 'Unknown' cell cycle
#' state prediction. ccAFv2 cell cycle state predictions and probabilities for
#' each cell in the object will be stored in the object .obs after classification.
#'
#'
#' @param Seurat object that should have ccAFv2 cell cycle states predicted.
#' @return Seurat object with ccAFv2 calls and probabilities for each cell cycle state.
#' @export
PredictCellCycle = function(seurat1, cutoff=0.5, assay='SCT', species='human', gene_id='ensembl') {
    cat('Running ccAFv2:\n')
    # Load model and marker genes
    ccAFv2 = keras::load_model_hdf5(system.file('extdata', 'ccAFv2_model.h5', package='ccAFv2'))
    classes = read.csv(system.file('extdata', 'ccAFv2_classes.txt', package='ccAFv2'), header=FALSE)[,1]
    mgenes = read.csv(system.file('extdata', 'ccAFv2_genes.csv', package='ccAFv2'), header=TRUE, row.names=1)[,paste0(species,'_',gene_id)]
    
    # Subset data marker genes to marker genes included in classification
    sub_genes = intersect(row.names(seurat1),mgenes)
    seurat_subset = subset(seurat1, features = sub_genes)
    
    # Find missing genes and assign 0s to each cell
    cat(paste0('  Total possible marker genes for this classifier: ', length(mgenes),'\n'))
    missing_genes = setdiff(mgenes, rownames(seurat_subset[[assay]]@data))
    cat(paste0('    Marker genes present in this dataset: ', nrow(seurat_subset[[assay]]@data),'\n'))
    cat(paste0('    Missing marker genes in this dataset: ', length(missing_genes),'\n'))
    ## Add ERROR later to warn if not enough marker genes ##
    input_mat = seurat_subset[[assay]]@data
    input_mat_scaled = t(scale(t(as.matrix(input_mat))))
    tmp = matrix(min(input_mat_scaled,na.rm=T),nrow=length(missing_genes), ncol=ncol(seurat1))
    rownames(tmp) = missing_genes
    colnames(tmp) = colnames(seurat_subset)
    input_mat_scaled_add_missing_genes = rbind(input_mat_scaled, tmp)[mgenes,]
    cat(paste0('  Predicting cell cycle state probabilities...\n'))
    predictions1 = predict(ccAFv2, t(input_mat_scaled_add_missing_genes))
    colnames(predictions1) = classes
    rownames(predictions1) = colnames(seurat1)
    df1 = data.frame(predictions1)
    cat(paste0('  Choosing cell cycle state...\n'))
    CellCycleState = data.frame(colnames(predictions1)[apply(predictions1,1,which.max)], row.names = rownames(predictions1))
    colnames(CellCycleState) = 'ccAFv2'
    df1[,'ccAFv2'] = CellCycleState$ccAFv2
    df1[which(apply(predictions1,1,max)<cutoff),'ccAFv2'] = 'Unknown'
    cat('  Adding probabilitities and predictions to metadata\n')
    seurat1 = AddMetaData(object = seurat1, metadata = df1)
    cat('Done\n')
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
    dp1 = DimPlot(seurat1, group.by='ccAFv2', cols = c('G1' = '#f37f73', 'G1/other' = '#9aca3c', 'G2/M' = '#3db270', 'Late G1' = '#1fb1a9','M/Early G1' = '#6d90ca', 'Neural G0' = '#d9a428', 'S' = '#8571b2', 'S/G2' = '#db7092'), ...)
    return(dp1)
}
