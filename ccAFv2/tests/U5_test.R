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
setwd('/files/U5_test')

# Download data file
download.file(url="https://zenodo.org/records/10961633/files/U5_normalized_ensembl.rds?download=1", destfile="U5_normalized_ensembl.rds", mode="wb")

# Load up sample data
seurat_obj = readRDS('U5_normalized_ensembl.rds')

# Run ccAFv2
seurat_obj = PredictCellCycle(seurat_obj,
                 threshold=0,
                 include_g0 = TRUE,
                 do_sctransform=TRUE,
                 assay='SCT',
                 layer = 'scale.data',
                 species='human',
                 gene_id='ensembl',
                 spatial = FALSE)

# Check to ensure 'ccAFv2' column exists
colnames(seurat_obj@meta.data)

# Show the first few rows of meta.data to see the likelihood estimates for each
# state, and ccAFv2 predictions in the 'ccAFv2' column
head(seurat_obj@meta.data)

# Number of cells in each cell cycle state
table(seurat_obj@meta.data$ccAFv2)

# Plot the UMAP with ccAFv2 predicitions
png('U5_hNSCs.png', width = 400, height = 400)
DimPlot.ccAFv2(seurat_obj)
dev.off()

