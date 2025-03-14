# ccAFv2: Cell cycle classifier for R and Seurat
[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)

This repository is for the R package for the cell cycle classifier ccAFv2. The input for the ccAFv2 classifier is single cell, nuclei, or spatial RNA-seq data.  The features of this classifier are that it classifies six cell cycle states (G1, Late G1, S, S/G2, G2/M, and M/Early G1) and a quiescent-like G0 state, and it incorporates a tunable parameter to filter out less certain classifications. This package is implemented in R for use in [Seurat](https://satijalab.org/seurat/) analysis workflows. We provide examples of installing and running ccAFv2 on Seurat objects (both sc/snRNA-seq and ST-RNA-seq), plot and use results, and regress out cell cycle effects if desired.

Installation descriptions are given below, and the API, tutorials, and more details please take a look at our github.io site:  https://plaisier-lab.github.io/ccafv2_R/

## Table of Contents

- [Install](#install)
    - [Requirements](#requirements)
        - [Dependencies](#dependencies)
        - [Dockerfile](#dockerfile)
        - [Docker images](#docker-images)
        - [Seurat 4.X version](#seurat-4.X-version)
        - [Seurat 5.X version](#seurat-5.X-version)
    - [Installing ccAFv2](#installing-ccafv2)
- [Tutorials](https://plaisier-lab.github.io/ccafv2_R/)
    - [Classifying single cell or nuclei RNA-seq](https://plaisier-lab.github.io/ccafv2_R/src/single.html)
    - [Classifying spatial RNA-seq](https://plaisier-lab.github.io/ccafv2_R/src/spatial.html)
    - [Cell cycle regression](https://plaisier-lab.github.io/ccafv2_R/src/regress.html)
    - [Choosing threshold](https://plaisier-lab.github.io/ccafv2_R/src/Choosing_Threshold.html)
- [Maintainers](#maintainers)
- [Contributing](#contributing)
- [Citation](#citation)

## Install

### Requirements

> It is strongly suggested that users utilize the Docker images we provide on DockerHub, as they contain all the dependencies needed to run ccAFv2.
However, for buliding from source note that ccAFv2 uses tensorflow 2.12. Reticulate will by default install python 3.12 wich only supports tensorflow 2.16. 
```r
reticulate::install_python(version = '3.9')
library(tensorflow)
install_tensorflow(
    method               = "auto", 
    version              = "2.12.1", 
    envname              = "retic", 
    conda_python_version = "3.9", 
    extra_packages       = c("matplotlib", "numpy", "pandas", "scikit-learn")
)
```

#### Dependencies

These dependencies must be met to run ccAFv2 classification:

- System dependencies:
    - Tensorflow/Keras version 2 must be installed - [install link](https://www.tensorflow.org/install)
- R dependencies:
    - Seurat 4 or 5 - [install link](https://satijalab.org/seurat/articles/install_v5)
    - keras - [install link](https://cran.r-project.org/web/packages/keras/vignettes/)

#### Dockerfile

We provide a Dockerfile that installs all the dependencies needed to run ccAFv2 so that it can be merged with a user's Dockerfile.

- [Dockerfile](https://github.com/plaisier-lab/ccafv2_R/blob/main/Dockerfile)

#### Docker images

We also provide fully compiled Docker images on DockerHub for ease of use. We provide two versions that provide users options as to which version of Seurat they use (4.X and 5.X).

##### Seurat 4.X version

Link to DockerHub image:  [cplaisier/ccafv2_seurat4](https://hub.docker.com/r/cplaisier/ccafv2_seurat4)

Command to pull the image down:

```
docker pull cplaisier/ccafv2_seurat4
```

Command to run the docker image. Note that the <replace with the location to your files> should be replaced with the path to your files that you want to be mounted onto the docker instance. The files can then be found in /files on the instance and locally on your computer in the path specified.

```sh
docker run -it -v '<replace with the location for your files>:/files' cplaisier/ccafv2_seurat4
```

##### Seurat 5.X version

Link to DockerHub image:  [cplaisier/ccafv2_seurat5](https://hub.docker.com/r/cplaisier/ccafv2_seurat5)

Command to pull the image down:

```
docker pull cplaisier/ccafv2_seurat5
```
Command to run the docker image. Note that the <replace with the location to your files> should be replaced with the path to your files that you want to be mounted onto the docker instance. The files can then be found in /files on the instance and locally on your computer in the path specified.

```sh
docker run -it -v '<replace with the location for your files>:/files' cplaisier/ccafv2_seurat5
```

### Installing ccAFv2

> **NOTE**: The Docker images already have ccAFv2 installed, so both of these commands are unnecessary if you use them.

Once the dependencies are met, ccAFv2 can be installed in R using the [remotes](https://cran.r-project.org/web/packages/remotes/index.html) package, which must be installed first. The devtools package can be installed using the command:

```r
install.packages('remotes')
```

Once the devtools package is installed, it can then be used to install the ccAFv2 R package from this GitHub repository using the following command:

```r
remotes::install_github('plaisier-lab/ccafv2_R/ccAFv2')
```

## Maintainers

For issues or comments, please contact:  [Chris Plaisier](mailto:plaisier@asu.edu)

For other great packages from the Plaisier Lab, please check here: [@plaisier-lab](https://github.com/plaisier-lab).

## Contributing

Feel free to dive in! [Open an issue](https://github.com/plaisier-lab/ccAFv2_R/issues/new) or submit PRs.

## Citation

1. **Citation for ccAFv2 (version 2)**:

[Classifying cell cycle states and a quiescent-like G0 state using single-cell transcriptomics](https://doi.org/10.1101/2024.04.16.589816) Samantha A. O’Connor, Leonor Garcia, Anoop P. Patel, Benjamin B. Bartelle, Jean-Philippe Hugnot, Patrick J. Paddison, Christopher L. Plaisier. bioRxiv [Preprint]. 2024 Apr 20:2024.04.16.589816. doi: [10.1101/2024.04.16.589816](https://doi.org/10.1101/2024.04.16.589816). PMID: [38659838](https://pubmed.ncbi.nlm.nih.gov/38659838/) 

2. **Citation for ccAF (version 1)**:

[Neural G0: a quiescent-like state found in neuroepithelial-derived cells and glioma.](https://www.embopress.org/doi/full/10.15252/msb.20209522) Samantha A. O'Connor, Heather M. Feldman, Chad M. Toledo, Sonali Arora, Pia Hoellerbauer, Philip Corrin, Lucas Carter, Megan Kufeld, Hamid Bolouri, Ryan Basom, Jeffrey Delrow, Jose L. McFaline-Figueroa, Cole Trapnell, Steven M. Pollard, Anoop Patel, Patrick J. Paddison, Christopher L. Plaisier.  Mol Syst Biol. 2021 Jun;17(6):e9522. doi: [10.15252/msb.20209522](https://doi.org/10.15252/msb.20209522). PMID: [34101353](https://pubmed.ncbi.nlm.nih.gov/34101353/)

