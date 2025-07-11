---
layout: default
title: Installation
nav_order: 3
---
# Requirements and dependencies

The ccAFv2 package is designed to work alongside [Seurat](https://satijalab.org/seurat), and makes predictions from and stores results into Seurat objects. It is suggested that the latest version of Seurat be installed. The classifier has been tested on Seurat version 5.3.1.

Installation of ccAFv2 from github is facilitated by the [remotes](https://cran.r-project.org/web/packages/remotes/index.html) package. Finally, for plotting we require installation of the [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) package.

```r
install.packages(c('Seurat', 'remotes', 'ggplot2'))
```

#  Installation

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


