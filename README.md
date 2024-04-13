# ccAFv2: Cell cycle classifier for R and Seurat
[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)

This repository is for the R package for the cell cycle classifier ccAFv2. The input for the ccAFv2 classifier is single cell or nuclei, or spatial RNA-seq data.  The features of this classifier are that it classifies six cell cycle states (G1, Late G1, S, S/G2, G2/M, and M/Early G1) and a quiescent-like G0 state, and it incorporates a tunable parameter to filter out less certain classifications. This package is implemented in R so that it can be used in [Seurat](https://satijalab.org/seurat/) analysis workflows. We provide examples of how to install, run ccAFv2 on seurat objects (both sc/snRNA-seq and ST-RNA-seq), plot and use results, and regress out cell cycle effects if that is desired.

## Table of Contents

- [Install](#install)
    - [Requirements](#requirements)
        - [Dependencies](#dependencies)
        - [Dockerfile](#dockerfile)
        - [Docker images](#docker-images)
        - [Seurat 4.X version](#seurat-4.X-version)
        - [Seurat 5.X version](#seurat-5.X-version)
    - [Installing ccAFv2](#installing-ccafv2)
- [Classifying single cell or nuclei RNA-seq](#classifying-single-cell-or-nuclei-rna-seq)
    - [Input for classification](#input-for-classification)
    - [Test data](#test-data)
	- [Cell cycle classification](#cell-cycle-classification)
    - [Plotting cell cycle states](#plotting-cell-cycle-states)
        - [Plotting a UMAP with cell cycle states](#plotting-a-umap-with-cell-cycle-states)
        - [Plotting the impact of varying likelihood thresholds](#plotting-the-impact-of-varying-likelihood-thresholds)
	- [Cell cycle regression](#cell-cycle-regression)
- [Classifying spatial RNA-seq](#classifying-spatial-rna-seq)
    - [Test spatial data](#test-spatial-data)
	- [Spatial cell cycle classification](#spatial-cell-cycle-classification)
    - [Plotting cell cycle states onto images](#plotting-cell-cycle-states-onto-images)
- [Maintainers](#maintainers)
- [Contributing](#contributing)
- [Citation](#citation)

## Install

### Requirements

> It is strongly suggested that users utilize the Docker images we provide on DockerHub, as they contain all the dependencies needed to run ccAFv2.

#### Dependencies

These dependencies must be met to run ccAFv2 classification:

- System dependencies:
    - Tensorflow/Keras version 2 must be installed - [install link](https://www.tensorflow.org/install)
- R dependencies:
    - Seurat 4 or 5 - [install link](https://satijalab.org/seurat/articles/install_v5)
    - keras - [install link](https://cran.r-project.org/web/packages/keras/vignettes/)

#### Dockerfile

We provide a Dockerfile that install all the dependencies needed to run ccAFv2 so that it can be merged with a users Dockerfile.

- [Dockerfile](https://github.com/plaisier-lab/ccafv2_R/blob/main/Dockerfile)

#### Docker images

We also provide fully compiled Docker images on DockerHub for ease of use. We provide two versions that provide users options as to which version of Seurat they are using (4.X and 5.X).

##### Seurat 4.X version

Link to DockerHub image:  [cplaisier/ccafv2_seurat4](https://hub.docker.com/r/cplaisier/ccafv2_seurat4)

Command to pull the image down:

```
docker pull cplaisier/ccafv2_seurat4
```

Command to run the docker image. Note that the <replace with the location to your files> should be replaced with the path to your files that you want mounted onto the docker instance. The files can then be found in /files on the instance, and locally on your computer in the path specified.

```sh
docker run -it -v '<replace with the location for your files>:/files' cplaisier/ccafv2_seurat4
```

##### Seurat 5.X version

Link to DockerHub image:  [cplaisier/ccafv2_seurat5](https://hub.docker.com/r/cplaisier/ccafv2_seurat5)

Command to pull the image down:

```
docker pull cplaisier/ccafv2_seurat5
```
Command to run the docker image. Note that the <replace with the location to your files> should be replaced with the path to your files that you want mounted onto the docker instance. The files can then be found in /files on the instance, and locally on your computer in the path specified.

```sh
docker run -it -v '<replace with the location for your files>:/files' cplaisier/ccafv2_seurat5
```

### Installing ccAFv2

> **NOTE**: The Docker images already have ccAFv2 installed, and this both of these commands are unnecessary if you use the Docker images.

Once the depnencies are met ccAFv2 can be installed in R using the [devtools](https://cran.r-project.org/web/packages/devtools/readme/README.html) package, which muse be installed first. The devtools packace can be installed using the command:

```r
install.packages('devtools')
```

Once the devtools package is installed it can then be used to install the ccAFv2 R package from this github repository using the following command:

```r
devtools::install_github('plaisier-lab/ccafv2_R/ccAFv2')
```

## Classifying single cell or nuclei RNA-seq

### Input for classification

It is expected that the input for the ccAFv2 classifier will be a Seurat object that has been thorougly quality controlled. We provide an example of our quality control pipeline in can be found [here](https://github.com/plaisier-lab/ccAFv2/blob/main/scripts/02_scQC_2024.R). Is is preferred that the data in the Seurat object be SCTransformed, however, the standard approach for normalization only applies to the highly variable genes. This can exclude genes needed for the
accurate classification of the cell cycle. For this reason the ccAFv2 PredictCellCycle function used to classify cell cycle states runs the SCTransform function again parameterized so that it will retain all genes captured in the dataset.

### Test data

The U5 human neural stem cell (hNSC) dataset used to train the ccAFv2 is available for testing purposes here:
- [U5 hNSCs rds file](https://zenodo.org/records/10961633/files/U5_normalized_ensembl.rds?download=1)

Download this file and place it into the directory in which you wish to run the ccAFv2 tutorial below. This data has been QC'd and normalized using SCTransform following our best practices described above.

### Cell cycle classification

Classification is as easy as two lines that can be added to any Seurat workflow. First the library must be loaded and then the PredictCellCycle function is run:

```r
library(ccAFv2)
seurat_obj = readRDS('U5_normalized_ensembl.rds')
seurat_obj = PredictCellCycle(seurat_obj)
```
When the classifier is running it should look something like this:

```r
Running ccAFv2:
 Redoing SCTransform to ensure maximum overlap with classifier genes...
 Total possible marker genes for this classifier: 861
  Marker genes present in this dataset: 861
  Missing marker genes in this dataset: 0
 Predicting cell cycle state probabilities...
93/93 [==============================] - 1s 4ms/step
93/93 [==============================] - 1s 4ms/step
 Choosing cell cycle state...
 Adding probabilitities and predictions to metadata
Done
```

It is important to look at how many marker genes were present in the dataset. We found that when less than 689 marker genes (or 80%) were found in the dataset that this led significantly less accurate predictions. Using the default 'do_sctransform' paramter setting of TRUE should yeild the largest possible overlap with the marker genes. And some of the later values for the timing and 93/93 may differ for your dataset, which is perfectly fine.

There are several options that can be passed to the PredictCellCycle function:
```r
PredictCellCycle(seurat_obj, 
                 cutoff=0.5, 
                 do_sctransform=TRUE,
                 assay='SCT',
                 species='human',
                 gene_id='ensembl',
                 spatial = FALSE) 
```
- **seurat_obj**: a seurat object must be supplied to classify, no default
- **cutoff**: the value used to threchold the likelihoods, default is 0.5
- **do_sctransform**: whether to do SCTransform before classifying, default is TRUE
- **assay**: which seurat_obj assay to use for classification, helpful if data is prenormalized, default is 'SCT'
- **species**: from which species did the samples originate, either 'human' or 'mouse', defaults to 'human'
- **gene_id**: what type of gene ID is used, either 'ensembl' or 'symbol', defaults to 'ensembl'
- **spatial**: whether the data is spatial, defaults to FALSE


### Cell cycle classification results

The results of the cell cycle classification is stored in the seurat object metadata. The likelihoods for each cell cycle state can be found with the labels of each cell cycle state ('Neural.G0', 'G1', 'Late.G1', 'S', 'S.G2', 'G2.M', and 'M.Early.G1') and the classification for each cell can be found int the 'ccAFv2'. Here are the first 10 rows of the U5-hNSC predictions:

```r
head(seurat_obj@meta.data)
```

Which returns the following:

```
                 orig.ident       nCount_RNA nFeature_RNA percent.mito
AAACATACTAACCG-1 SeuratProject    3147       1466         0.05274865
AAACATTGAGTTCG-1 SeuratProject    7621       2654         0.02611206
AAACATTGCACTGA-1 SeuratProject    7297       2616         0.02686035
AAACATTGCTCAGA-1 SeuratProject    3426       1568         0.02597782
AAACATTGGTTTCT-1 SeuratProject    2384       1269         0.01971477
AAACCGTGTAACGC-1 SeuratProject   11846       3463         0.02439642
AAACGCACCTTCTA-1 SeuratProject    3021       1249         0.01886792
AAACGCTGGTATGC-1 SeuratProject    5937       2127         0.02846555
AAACGCTGTGCTGA-1 SeuratProject    8361       2793         0.03564167
AAACGGCTGTCTAG-1 SeuratProject    4591       1605         0.04530603
                 nCount_SCT nFeature_SCT G1           G2.M         Late.G1
AAACATACTAACCG-1 5436       1541         4.102312e-02 5.279975e-05 2.521812e-04
AAACATTGAGTTCG-1 6562       2650         1.091572e-06 2.126902e-08 1.087603e-07
AAACATTGCACTGA-1 6533       2616         4.089213e-06 1.517231e-04 6.863667e-07
AAACATTGCTCAGA-1 5416       1608         9.184604e-03 2.376643e-01 1.454412e-02
AAACATTGGTTTCT-1 5437       1480         6.363088e-01 1.654975e-04 9.073097e-02
AAACCGTGTAACGC-1 6860       3082         5.167215e-15 9.999999e-01 1.610888e-19
AAACGCACCTTCTA-1 5537       1308         4.707390e-06 3.709571e-08 6.463646e-05
AAACGCTGGTATGC-1 5954       2127         2.741130e-01 2.113478e-05 7.243306e-01
AAACGCTGTGCTGA-1 6666       2780         1.320201e-04 2.646190e-07 1.700663e-06
AAACGGCTGTCTAG-1 5568       1608         3.885373e-01 3.646265e-03 2.107647e-02
                 M.Early.G1   Neural.G0    S            S.G2         ccAFv2
AAACATACTAACCG-1 4.220059e-05 9.585268e-01 5.953405e-05 4.332590e-05 Neural G0
AAACATTGAGTTCG-1 8.592179e-10 9.998107e-01 1.868801e-04 1.347080e-06 Neural G0
AAACATTGCACTGA-1 4.317539e-09 1.362584e-06 6.708771e-07 9.998415e-01 S/G2
AAACATTGCTCAGA-1 6.309561e-01 5.861283e-02 4.330526e-03 4.470748e-02 M/Early G1
AAACATTGGTTTCT-1 2.112151e-05 1.050980e-03 1.899571e-01 8.176555e-02 G1
AAACCGTGTAACGC-1 1.919447e-18 3.083377e-16 6.074685e-23 1.652127e-15 G2/M
AAACGCACCTTCTA-1 4.957183e-07 9.401034e-07 9.998928e-01 3.640788e-05 S
AAACGCTGGTATGC-1 7.537001e-06 3.954049e-04 9.016516e-04 2.307067e-04 Late G1
AAACGCTGTGCTGA-1 1.102674e-08 9.998627e-01 2.419535e-06 1.108445e-06 Neural G0
AAACGGCTGTCTAG-1 7.390769e-04 5.405837e-01 5.436572e-03 3.998061e-02 Neural G
```

### Plotting cell cycle states

We provide plotting functions that colorize the cell cycle states in the way used in our manuscripts. We strongly suggest using these functions when plotting if possible.

#### Plotting a UMAP with cell cycle states

Plotting cells using ther first two dimensions from a dimensionality reduction method (e.g., PCA, tSNE, or UMAP) is a common way to represent single cell or nuclei RNA-seq data. We have an overloaded DimPlot function that colorizes the cells based on their called cell cycle state. The function accepts all the parameters that DimPlot can accept, except for group.by and cols. Here is how the plotting function should be run:

```r
pdf('ccAFv2_DimPlot.pdf')
DimPlot.ccAFv2(seurat_obj)
dev.off()
```
Below is a DimPlot for U5 hNSCs colorized using the cell cycle states. The expected flow of the cell cycle states can be seen in the UMAP.

![UMAP DimPlot colorized with ccAFv2 cell cycle states](https://github.com/plaisier-lab/ccAFv2_R/blob/main/figs/ccAFv2_DimPlot.png?raw=true)

#### Plotting the impact of varying likelihood thresholds

For certain datasets it may be necessary to adjust the likelihood threshold. We provide a plot that can make this process easier. The ThresholdPlot plots the relative proportions of cell cycle states across a range or likelihood thresholds from 0 to 0.9 by intervals of 0.1.

```r
pdf('ccAFv2_ThresholdPlot.pdf')
ThresholdPlot(seurat_obj)
dev.off()
```

Likelihood thersholds are on the x-axis and relative proportions of cell cycle states are on the y-axis. As can be seen as the likelihood thresholds increase the number of 'Unknown' cells increases.

![Bar plot colorized with ccAFv2 cell cycle states acorss varios thresholds](https://github.com/plaisier-lab/ccAFv2_R/blob/main/figs/ccAFv2_ThresholdPlot.png?raw=true)

### Cell cycle regression

The cell cycle imposes a strong biological signal on cell expression patterns. Thus it has become a common practice to regress the cell cycle expression out of cells expression, and then use the residual variance for further studies. We provide functionality to do this using the ccAFv2 marker genes.

First, we collect expression module scores for the cell cycle states Late G1, S, S/G2, G2/M, and M/Early G1.

```r
seurat_obj = PrepareForCellCycleRegression(seurat_obj)
```

Then we regress these signatures out of the expression data:

```r
seurat_obj = SCTransform(seurat_obj, vars.to.regress = c("Late.G1_exprs1", "S_exprs2", "S.G2_exprs3", "G2.M_exprs4", "M.Early.G1_exprs5"))
```

And finally to plot the effect of regressing out the cell cycle on the UMAP:

```r
seurat_obj = RunPCA(seurat_obj)
seurat_obj = RunUMAP(seurat_obj, dims=1:10)
pdf('ccAFv2_DimPlot_regressed.pdf')
DimPlot.ccAFv2(seurat_obj)
dev.off()
```

Removing the cell cycle from the U5 hNSCs leads to a random distribution, because the cell cycle is the primary biological signal in these *in vitro* grown cell line.

![UMAP DimPlot colorized with ccAFv2 cell cycle states after regressing out cell cycle](https://github.com/plaisier-lab/ccAFv2_R/blob/main/figs/ccAFv2_DimPlot_regressed.png?raw=true)


## Classifying spatial RNA-seq

### Test spatial data

A slice of a human fetus 4 weeks post conception from [Zeng et al., 2023](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/geo/query/acc.cgi?acc=GSE155121) is available for testing purposes here:
- [Zeng et al., 2023 - 4 weeks post conception human fetus - slice 1](https://zenodo.org/records/10961633/files/GSM6736780_Spatial_10x_PCW4_20220122_slice1.rds?download=1)

Download this file and place it into the directory in which you wish to run the ccAFv2 spatial tutorial below. This data has been QC'd and normalized using SCTransform following our best practices described above.

### Cell cycle classification

Classification is as easy as two lines that can be added to any Seurat workflow. First the library must be loaded and then the PredictCellCycle function is run:

```r
library(ccAFv2)
spatial_obj = readRDS('GSM6736780_Spatial_10x_PCW4_20220122_slice1.rds')
spatial_obj = PredictCellCycle(spatial_obj, species='human', gene_id='symbol', spatial=TRUE)
```
When the classifier is running it should look something like this:

For details about expected output please see classifying cells and nuclei above.

### Plotting cell cycle states onto images

Plotting cell cycle states onto the images taken of the tissue slices before spatial RNA-seq is a common way to represent spatial RNA-seq data. We have an overloaded SpatialDimPlot function that colorizes the cells based on their called cell cycle state. The function accepts all the parameters that SpatialDimPlot can accept, except for group.by and cols. Here is how the plotting function should be run:

```r
pdf('ccAFv2_SpatialDimPlot_slice1.pdf')
SpatialDimPlot.ccAFv2(spatial_obj) + theme(legend.position = "right")
dev.off()
```
Below is a SpatialDimPlot for slice 1 of a human fetus at 4 weeks post-conception colorized using the cell cycle states.

![UMAP SpatialDimPlot colorized with ccAFv2 cell cycle states](https://github.com/plaisier-lab/ccAFv2_R/blob/main/figs/ccAFv2_SpatialDimPlot_slice1.png?raw=true)

## Maintainers
For issues or comments please contact:  [Chris Plaisier](mailto:plaisier@asu.edu)
And for other great packages from the Plaisier Lab please check here:  [@plaisier-lab](https://github.com/plaisier-lab).

## Contributing

Feel free to dive in! [Open an issue](https://github.com/plaisier-lab/ccAFv2_R/issues/new) or submit PRs.

## Citation

> - **Citation for ccAFv2 (version 2)**:
>
>   *Citation for ccAFv2 coming soon!*

- **Citation for ccAF (version 1)**:
[Neural G0: a quiescent-like state found in neuroepithelial-derived cells and glioma.](https://doi.org/10.1101/446344) Samantha A. O'Connor, Heather M. Feldman, Chad M. Toledo, Sonali Arora, Pia Hoellerbauer, Philip Corrin, Lucas Carter, Megan Kufeld, Hamid Bolouri, Ryan Basom, Jeffrey Delrow, Jose L. McFaline-Figueroa, Cole Trapnell, Steven M. Pollard, Anoop Patel, Patrick J. Paddison, Christopher L. Plaisier. bioRxiv 446344; doi: [https://doi.org/10.1101/446344](https://doi.org/10.1101/446344)
