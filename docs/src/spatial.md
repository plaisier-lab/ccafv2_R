---

layout: default
title: Classifying Spatial RNA-seq
nav_order: 10 

---
# Input For Classification

The ccAFv2 classifier requires a thoroughly quality-controlled Seurat
object. SCTransform is preferred, standard normalization only applies to the
highly variable genes. This can exclude genes needed for the accurate
classification of the cell cycle. To address this, the PredictCellCycle
function re-runs SCTransform to retain all genes in the dataset

# Test Spatial Data

A slice of a mouse feuts from 10x Genomics who obtained FFPE Mouse Embryo tissue blocks from Avaden Biosciences. Required the download of the [images](https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_11mm_FFPE_Mouse_Embryo/CytAssist_11mm_FFPE_Mouse_Embryo_spatial.tar.gz) and [h5 file](https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_11mm_FFPE_Mouse_Embryo/CytAssist_11mm_FFPE_Mouse_Embryo_filtered_feature_bc_matrix.h5). This data has been QC'd and normalized using SCTransform, following our
best practices described above.

## Load libraries
```r
library(ccAFv2)
library(Seurat)
```

## Load Data
```r
image1 = Read10X_Image(image.dir = 'spatial')
spatial_obj = Load10X_Spatial(data.dir = 'filtered_feature_bc_matrix', image=image1)
```

## Scale and Normalize Data
```r
spatial_obj = SCTransform(spatial_obj, assay='Spatial', verbose=FALSE)
```

## Predict Cell Cycle States
Note that this particular samples used gene symbols, which is different from our default of Ensembl IDs. Which is why that is specified in the functional call.
```r
spatial_obj = PredictCellCycle(spatial_obj,
                               include_g0 = TRUE, 
                               species='mouse', 
                               gene_id='symbol', 
                               spatial=TRUE)
```

# Plotting Cell Cycle States Onto Images

To represent spatial RNA-seq data, cell cycle states are often mapped
onto tissue slice images taken before sequencing. We use an enhanced
**SpatialDimPlot** function that colors cells by their cell cycle state.
This function supports all **SpatialDimPlot** parameters except
**group.by** and **cols**.

```r
pdf('spatial_ccAFv2.pdf')
SpatialPlot.ccAFv2(spatial_obj)
dev.off()
```
## Combined G0/G1
SpatialDimPlot for slice 1 of a mouse fetus colorized using the cell cycle states that combine G0 and G1 states.

![ccAFv2_SpatialDimPlot_G0_G1]({{ site.baseurl }}/images/spatial_dim_plot.png)

## Separate G0 and G1

SpatialDimPlot for slice 1 of a mouse fetus colorized using the cell cycle states that combine G0 and G1 states.

![ccAFv2_SpatialDimPlot_G0]({{ site.baseurl }}/images/spatial_dim_plot_WG0.png)

## References 
Zeng B, Liu Z, Lu Y, Zhong S, Qin S, Huang L, Zeng Y, Li Z, Dong H, Shi Y, et al. 2023. *The single-cell and spatial transcriptional landscape of human gastrulation and early brain development.* Cell Stem Cell 30: 851–866.e7. [PubMed](https://pubmed.ncbi.nlm.nih.gov/37192616/) 

## Load Data
```r
spatial_obj = readRDS('GSM6736780_Spatial_10x_PCW4_20220122_slice1.rds')
```

## Predict Cell Cycle States
Note that this particular samples used gene symbols, which is different from our default of Ensembl IDs. Which is why that is specified in the functional call.
```r
spatial_obj = PredictCellCycle(spatial_obj, species='human', gene_id='symbol', spatial=TRUE)
```

# Plotting Cell Cycle States Onto Images

To represent spatial RNA-seq data, cell cycle states are often mapped
onto tissue slice images taken before sequencing. We use an enhanced
**SpatialDimPlot** function that colors cells by their cell cycle state.
This function supports all **SpatialDimPlot** parameters except
**group.by** and **cols**.

```r
library(ggplot2)
SpatialDimPlot.ccAFv2(spatial_obj) + theme(legend.position = "right")

```

SpatialDimPlot for slice 1 of a human fetus at 4 weeks post-conception
colorized using the cell cycle states.

![ccAFv2_SpatialDimPlot_slice1]({{ site.baseurl }}/images/ccAFv2_SpatialDimPlot_slice1.png)

## References 
Zeng B, Liu Z, Lu Y, Zhong S, Qin S, Huang L, Zeng Y, Li Z, Dong H, Shi Y, et al. 2023. *The single-cell and spatial transcriptional landscape of human gastrulation and early brain development.* Cell Stem Cell 30: 851–866.e7. [PubMed](https://pubmed.ncbi.nlm.nih.gov/37192616/) 
