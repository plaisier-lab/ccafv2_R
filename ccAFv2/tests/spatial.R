##########################################################
## ccAFv2:  Test on U5 hNSCs                            ##
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
## @Author:  Chris Plaisier                             ##
## @License:  GNU GPLv3                                 ##
##                                                      ##
## If this program is used in your analysis please      ##
## mention who built it. Thanks. :-)                    ##
##########################################################

# Get docker to run analyses
# docker pull satijalab/seurat ## R V4.6.1, Seurat V6.0.0, SeuratObject V5.4.0, 24.04.5 LTS (Noble Numbat), base image rocker/tidyverse:4.6.1
# docker run -it -v '/media/old_home/home/cplaisier/ccAFv2_test:/files' satijalab/seurat

# Install required packages
install.packages(c('remotes','ggplot2', 'circular'))
remotes::install_github('plaisier-lab/ccafv2_R/ccAFv2')

# Load libraries
library(ccAFv2)
library(Seurat)

# Set working directory
setwd('/files/Mm_embryo')

# Load up data
image1 = Read10X_Image(image.dir = 'spatial')
spatial_obj = Load10X_Spatial(data.dir = 'filtered_feature_bc_matrix', image=image1)

# Scale and normalize data
spatial_obj = SCTransform(spatial_obj, assay='Spatial', verbose=FALSE)

# Run ccAFv2
spatial_obj = PredictCellCycle(spatial_obj,
                               include_g0 = F, 
                               species='mouse', 
                               gene_id='symbol', 
                               spatial=TRUE)

# Plot SpatialDimPlot
png('spatial_dim_plot.png', height = 800, width = 800)
p = SpatialDimPlot.ccAFv2(spatial_obj)
print(p)
dev.off()

# Run ccAFv2
spatial_obj = PredictCellCycle(spatial_obj,
                               include_g0 = TRUE, 
                               species='mouse', 
                               gene_id='symbol', 
                               spatial=TRUE)

# Plot SpatialDimPlot
png('spatial_dim_plot_wG0.png', height = 800, width = 800) 
p = SpatialDimPlot.ccAFv2(spatial_obj)
print(p)
dev.off()

