##########################################################
## ccAFv2: project_cell_cycle.R                         ##
##  ______     ______     __  __                        ##
## /\  __ \   /\  ___\   /\ \/\ \                       ##
## \ \  __ \  \ \___  \  \ \ \_\ \                      ##
##  \ \_\ \_\  \/\_____\  \ \_____\                     ##
##   \/_/\/_/   \/_____/   \/_____/                     ##
## @Developed by: Plaisier Lab                          ##
##   (https://plaisierlab.engineering.asu.edu/)         ##
##   Arizona State University                           ##
##   242 ISTB1, 550 E Orange St                         ##
##   Tempe, AZ  85281                                   ##
## @Author:  Chris Plaisier, Samantha O'Connor          ##
## @License:  GNU GPLv3                                 ##
##                                                      ##
## If this program is used in your analysis please      ##
## mention who built it. Thanks. :-)                    ##
##########################################################

## This is derived from and credit is given to https://github.com/hansenlab/tricycle,
#  major modification is integration with Seurat and the use of Seurat SCT
#  transformed and normalized data.

## TODO - Gene ID conversions, species conversions

#' @importFrom methods is
#' @importClassesFrom AnnotationDbi AnnotationDb
.GetAnnotationDB = function(AnnotationDB = NULL, species = c("mouse", "human")) {
    species = match.arg(species)
    if (!is.null(AnnotationDB)) {
        stopifnot("The input AnnotationDB is not an AnnotationDb object.", is(AnnotationDB, "AnnotationDb"))
    } else {
        if (species == "mouse") {
            stopifnot("No AnnotationDB input and org.Mm.eg.db package is not installed. Please install org.Mm.eg.db package." = nzchar(system.file(package = "org.Mm.eg.db")))
            AnnotationDb = org.Mm.eg.db::org.Mm.eg.db
            message("No AnnotationDb desginated. org.Mm.eg.db will be used to map Mouse ENSEMBL id to gene SYMBOL.")
        } else {
            stopifnot("No AnnotationDB input and org.Hs.eg.db package is not installed. Please install org.Mm.Hs.db package." = nzchar(system.file(package = "org.Hs.eg.db")))
            AnnotationDb = org.Hs.eg.db::org.Hs.eg.db
            message("No AnnotationDb desginated. org.Hs.eg.db will be used to map Human ENSEMBL id to gene SYMBOL.")
        }
    }
    return(AnnotationDb)
}

##' Get the expression data from a Seurat object
#' 
#' This function will collect the gene expression data from a Seurat object.
#'
#' @param seurat_obj: the Seurat object.
#' @param assay: which assay in the Seurat object should be exported, default is 'SCT'.
#' @param layer: which layer in the Seurat object assay should be exported, default is 'data'.
#' @param gname: list of gene ids as input, or use Seurat rownames if NULL, default is NULL.
#' @param gname.type: gene id type, either ENSEMBL or SYMBOL, no default.
#' @param species: whether using human or mouse, no default.
#' @param ntop: number of top cell cycle genes to include in the analysis, default is 500.
#' @param ncomponents: numer of principal components to detect, default is 20.
#' @param name: name of the PCA result in the seuratObject, default is 'cyclePCA'.
#' @param AnnotationDb: the annotation database, either org.Hs.eg.db or org.Mm.eg.db, default is NULL.
#' @return Seurat expression data.
#' @export
#' 
RunPCACCGenes = function(seurat_obj, 
                            assay='SCT',
                            layer='scale.data',
                            gname = NULL,
                            gname.type = c("ENSEMBL", "SYMBOL"),
                            species = c("mouse", "human"),
                            ntop = 500,
                            ncomponents = 20,
                            name = "cyclePCA",
                            AnnotationDb = NULL) {

  species = match.arg(species)
  gname.type = match.arg(gname.type)
  AnnotationDb = .GetAnnotationDB(AnnotationDb, species)
  if (is.null(gname)) {
    message("No gname input. Rownames of Seurat object will be used.")
    gname = rownames(seurat_obj)
  } else {
    stopifnot("gname does not match nrow Seruat object!." = nrow(seurat_obj) == length(gname))
  }

  cycle.anno = AnnotationDbi::select(AnnotationDb, keytype = "GOALL", keys = "GO:0007049", columns = gname.type)[, gname.type]
  gene.idx = which(gname %in% cycle.anno)
  if (length(gene.idx) < 100) {
    stop("Less than 100 Gene Ontology cell cycle genes found in your data. Check you data or gname input.")
  }
  message(paste0(length(gene.idx), " out of ", length(cycle.anno), " Gene Ontology cell cycle genes found in your data."))
  expr  = get_seurat_expr(seurat_obj, assay = assay, layer = layer)
  expr = expr[gene.idx, ]
  pca1 = runPCA(expr, ntop = ntop, ncomponents = ncomponents)
  return(pca1[,1:2])
}


##' Get the expression data from a Seurat object
#' 
#' This function will collect the gene expression data from a Seurat object.
#'
#' @param seurat_obj: the Seurat object.
#' @param assay: which assay in the Seurat object should be exported, default is 'SCT'.
#' @param layer: which layer in the Seurat object assay should be exported, default is 'data'.
#' @return Seurat expression data.
#' @export
#' 
GetSeuratExpr = function(seurat_obj, assay = 'SCT', layer = 'scale.data') {
  if (packageVersion('SeuratObject') >= '5.0.0') {
     # Seurat v5
     return(LayerData(
         object = seurat_obj,
         assay = assay,
         layer = layer
     ))
  } else {
     # Seurat v4
     return(GetAssayData(
         object = seurat_obj,
         assay = assay,
         slot = layer
     ))
  }
}


##' Funtion to project new cells into a reference embedding
#' 
#' Project cell cycle onto cells from an expression matrix.
#'
#' @param expr: an expression matrix.
#' @param ref.m: reference embedding on which to project the new cells.
#' @param center.pc1: center of PC1.
#' @param center.pc2: center of PC2.
#' @return Seurat expression data.
#' @export
#' 
ProjectCycleFromMatrix = function(expr, ref.m, center.pc1 = 0, center.pc2 = 0) {
  genes = intersect(rownames(expr), rownames(ref.m))
  if (length(genes) == 0) {
    stop('No overlapping genes between the Seurat assay and the reference CSV.')
  }
  message(paste0("    The number of projection genes found in the new data is ", length(genes), "."))

  ref.m = ref.m[genes, , drop = FALSE]
  expr  = expr[genes,  , drop = FALSE]
  proj = scale(t(as.matrix(expr)), center = TRUE, scale = FALSE) %*% ref.m
  rownames(proj) = colnames(expr)
  colnames(proj) = colnames(ref.m)
  attr(proj, "rotation") = ref.m

  proj[,1] = proj[,1] - center.pc1
  proj[,2] = proj[,2] - center.pc2
  theta = as.numeric(circular::coord2rad(proj[, seq_len(2)]))
  names(theta) = colnames(expr)

  list(projection = proj, position = theta)
}


##' Funtion to project new cells into a reference embedding and put it back into a Seurat object
#' 
#' Project cell cycle onto cells from an expression matrix and store
#' results in a Seurat object
#'
#' @param expr: an expression matrix.
#' @param ref_rds: name of the CSV file for the reference.
#' @param gene_col: which gene column should be used.
#' @param assay: which assay in the Seurat object should be exported, default is 'SCT'.
#' @param layer: which layer in the Seurat object assay should be exported, default is 'data'.
#' @param reduction_name: name of the reduction we generate using the projection.
#' @param reduction_key: string that should be prepended to reduction variable names.
#' @param position_col: column name for the storage of the theta angle for each cell in the metadata.
#' @param center.pc1: center of PC1.
#' @param center.pc2: center of PC2.
#' @return Seurat object with the projection data.
#' @export
#' 
ProjectCycleFromSeurat = function(
  seurat_obj,
  ref_rds,
  assay = 'SCT',
  layer = 'scale.data',
  gene_col = 1,
  reduction_name = 'cc_Projection',
  reduction_key = 'CYCLE_',
  position_col = 'thetaPos',
  center.pc1 = 0,
  center.pc2 = 0
) {
  stopifnot(inherits(seurat_obj, 'Seurat'))

  ref.m = ref_rds

  # Properly normalize and scale for projection
  expr = GetSeuratExpr(seurat_obj, assay = assay, layer = layer)

  res = ProjectCycleFromMatrix(
    expr = expr,
    ref.m = ref.m,
    center.pc1 = center.pc1,
    center.pc2 = center.pc2
  )
  
  if (!requireNamespace('SeuratObject', quietly = TRUE)) {
    stop('SeuratObject is required to store the reduction back into the object.')
  }

  seurat_obj[[reduction_name]] = SeuratObject::CreateDimReducObject(
    embeddings = as.matrix(res$projection),
    key = reduction_key,
    assay = assay
  )
  
  rotation = 0.225 * pi
  seurat_obj[[position_col]] = unlist(sapply(res$position, function(x) { (2*pi) - ((x - rotation) %% (2 * pi)) }))

  seurat_obj
}
