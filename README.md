# ccAFv2: Cell cycle classifier for R and Seurat
[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)
[![C](https://img.shields.io/badge/C-A8B9CC.svg?style=for-the-badge&logo=C&logoColor=black)](https://en.wikipedia.org/wiki/C_(programming_language))
[![R](https://img.shields.io/badge/R-276DC3.svg?style=for-the-badge&logo=R&logoColor=white)](https://cran.r-project.org/)
[![Docker](https://img.shields.io/badge/Docker-2496ED.svg?style=for-the-badge&logo=Docker&logoColor=white)](https://www.docker.com/)
[<img src="https://img.shields.io/badge/dockerhub-ccafv2-blue.svg?logo=Docker">](https://hub.docker.com/r/cplaisier/ccafv2)

This repository is for the R package for the cell cycle classifier ccAFv2. The input for the ccAFv2 classifier is single cell, nuclei, or spatial RNA-seq data.  The features of this classifier are that it classifies six cell cycle states (G1, Late G1, S, S/G2, G2/M, and M/Early G1) and a quiescent-like G0 state, and it incorporates a tunable parameter to filter out less certain classifications. This package is implemented in R for use in [Seurat](https://satijalab.org/seurat/) analysis workflows. We provide examples of installing and running ccAFv2 on Seurat objects (both sc/snRNA-seq and ST-RNA-seq), plot and use results, and regress out cell cycle effects if desired.

Installation descriptions are given below, and the API, tutorials, and more details please take a look at our github.io site:  https://plaisier-lab.github.io/ccafv2_R/

## Table of Contents

- [Requirements and dependencies](#requirements-and-dependencies)
- [Install](#install)
- [Dockerfile and image](#dockerfile-and-image)
- [Tutorials](https://plaisier-lab.github.io/ccafv2_R/)
    - [Classifying single cell or nuclei RNA-seq](https://plaisier-lab.github.io/ccafv2_R/src/single.html)
    - [Classifying spatial RNA-seq](https://plaisier-lab.github.io/ccafv2_R/src/spatial.html)
    - [Cell cycle regression](https://plaisier-lab.github.io/ccafv2_R/src/regress.html)
    - [Choosing threshold](https://plaisier-lab.github.io/ccafv2_R/src/Choosing_Threshold.html)
- [Maintainers](#maintainers)
- [Contributing](#contributing)
- [Citation](#citation)

## Requirements and dependencies

The ccAFv2 package is designed to work alongside [Seurat](https://satijalab.org/seurat), and makes predictions from and stores results into Seurat objects. It is suggested that the latest version of Seurat be installed. The classifier has been tested on Seurat version 5.3.1.

Installation of ccAFv2 from github is facilitated by the [remotes](https://cran.r-project.org/web/packages/remotes/index.html) package. Finally, for plotting we require installation of the [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) package.

```r
install.packages(c('Seurat', 'remotes', 'ggplot2'))
```

## Install

Installation of ccAFv2 in R is very simple and has very few requirements:

```r
remotes::install_github('plaisier-lab/ccafv2_R/ccAFv2')
```

## Dockerfile and image

For those who are interested we also provide a Dockerfile and image that include all depdencies:

- [Dockerfile](https://github.com/plaisier-lab/ccafv2_R/blob/main/Dockerfile)
- [DockerHub image:  cplaisier/ccafv2](https://hub.docker.com/r/cplaisier/ccafv2)

Command to pull the image down:

```sh
docker pull cplaisier/ccafv2
```

Command to run the docker image. Note that the <replace with the location to your files> should be replaced with the path to your files that you want to be mounted onto the docker instance. The files can then be found in /files on the instance and locally on your computer in the path specified.

```sh
docker run -it -v '<replace with the location for your files>:/files' cplaisier/ccafv2
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

