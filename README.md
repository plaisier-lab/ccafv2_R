# ccAFv2_R
[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)
This repository is for the R package for the cell cycle classifier ccAFv2. The input for the ccAFv2 classifier is single cell or nuclei, or spatial RNA-seq data.  The features of this classifier are that it classifies six cell cycle states (G1, Late G1, S, S/G2, G2/M, and M/Early G1) and a quiescent-like G0 state, and it incorporates a tunable parameter to filter out less certain classifications. This package is implemented in R so that it can be used in [Seurat](https://satijalab.org/seurat/) analysis workflows. We provide examples of how to install, run ccAFv2 on seurat objects (both sc/snRNA-seq and ST-RNA-seq), plot and use results, and regress out cell cycle effects if that is desired.

## Table of Contents

- [Install](#install)
- [Classifying single cell or nuclei RNA-seq](#classifying-single-cell-or-nuclei-rna-seq)
	- [Cell cycle classification](#cell-cycle-classification)
	- [Cell cycle likelihoods](#cell-cycle-classification)
    - [Plotting cell cycle states](#plotting-cell-cycle-states)
	- [Applying thresholds](#applying-thresholds)
	- [Cell cycle regression](#cell-cycle-regression)
- [Usage spatial RNA-seq](#usage)
	- [Cell cycle classification](#cell-cycle-classification)
	- [Cell cycle likelihoods](#cell-cycle-classification)
    - [Plotting cell cycle states](#plotting-cell-cycle-states)
	- [Applying thresholds](#applying-thresholds)
- [Maintainers](#maintainers)
- [Contributing](#contributing)

## Install

The installation of ccAv2 in R requires the [devtools](https://cran.r-project.org/web/packages/devtools/readme/README.html) package be install first. The devtools packace can be accomplished using the command:

```r
install.packages('devtools')
```

Once the devtools package is installed it can then be used to install the ccAFv2 R package using the following command:

```r
devtools::install_github('plaisier-lab/ccafv2_R/ccAFv2')
```

## Classifying single cell or nuclei RNA-seq

### Input for classification

It is expected that the input for the ccAFv2 classifier will be a Seurat object that has been thorougly quality controlled. We provide an example of our quality control pipeline in can be found [here](https://github.com/plaisier-lab/ccAFv2/blob/main/scripts/02_scQC_2024.R). Is is preferred that the data in the Seurat object be SCTransformed, however, the standard approach for normalization only applies to the highly variable genes. This can exclude genes needed for the
accurate classification of the cell cycle. For this reason the ccAFv2 PredictCellCycle function used to classify cell cycle states runs the SCTransform function again parameterized so that it will retain all genes captured in the dataset.

### Cell cycle classification

Classification is as easy as two lines that can be added to any Seurat workflow. First the library must be loaded and then the PredictCellCycle function is run:

```r
library(ccAFv2)
seurat_obj = PredictCellCycle(seurat_obj)

```
There are several options that can be passed to the PredictCellCycle function:
```
PredictCellCycle = function(seurat0, cutoff=0.5, do_sctransform=TRUE, assay='SCT', species='human', gene_id='ensembl', spatial = FALSE) 
```
- seurat0: seurat_obj, no default
- cutoff: the value used to threchold the likelihoods, default is 0.5
- do_sctransform: whether to do SCTransform before classifying, default is TRUE
- assay: which seurat_obj assay to use for classification, helpful if data is prenormalized, default is 'SCT'
- species: from which species did the samples originate, either 'human' or 'mouse', defaults to 'human'
- gene_id: what type of gene ID is used, either 'ensembl' or 'symbol', defaults to 'ensembl'
- spatial: whether the data is spatial, defaults to FALSE

### Cell cycle classification results



### Cell cycle regression

- [open-source-template](https://github.com/davidbgk/open-source-template/) - A README template to encourage open-source contributions.

## Maintainers

[@plaisier-lab](https://github.com/plaisier-lab).

## Contributing

Feel free to dive in! [Open an issue](https://github.com/plaisier-lab/ccAFv2_R/issues/new) or submit PRs.
