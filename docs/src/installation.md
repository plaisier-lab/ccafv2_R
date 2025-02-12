---
layout: default
title: Installation
nav_order: 3
---
#  Installation
It is strongly suggested that users utilize the Docker images we provide on DockerHub, as they contain all the dependencies needed to run ccAFv2
Installing ccAFv2
Contact the authors if you want to use an anconda virtual envrionment on your local machine to avoid dependency issues. 

NOTE: The Docker images already have ccAFv2 installed, so both of these commands are unnecessary if you use them.

Once the dependencies are met, ccAFv2 can be installed in R using the devtools package, which must be installed first. The devtools package can be installed using the command:
```
install.packages('devtools')
devtools::install_github('plaisier-lab/ccafv2_R/ccAFv2')
```
# System Dependencies
- Tensorflow/Keras version 2 must be installed - [install link](https://www.tensorflow.org/install)

R dependencies:
- Seurat 4 or 5 - [install link](https://satijalab.org/seurat/articles/install_v5)
- keras - [install link](https://cran.r-project.org/web/packages/keras/vignettes/)

# Dockerfile

Seurat 4.X version *Link to DockerHub image: [cplaisier/ccafv2_seurat4](https://hub.docker.com/r/cplaisier/ccafv2_seurat4)*
```
docker pull cplaisier/ccafv2_seurat4
```
Note that the should be replaced with the path to your files that you want to be mounted onto the docker instance. The files can then be found in /files on the instance and locally on your computer in the path specified.

```
docker run -it -v '<replace with the location for your files>:/files' cplaisier/ccafv2_seurat4
```
Seurat 5.X version *Link to DockerHub image: [cplaisier/ccafv2_seurat4](https://hub.docker.com/r/cplaisier/ccafv2_seurat5)*
```
docker pull cplaisier/ccafv2_seurat5
```
Command to run the docker image. Note that the should be replaced with the path to your files that you want to be mounted onto the docker instance. The files can then be found in /files on the instance and locally on your computer in the path specified.
```
docker run -it -v '<replace with the location for your files>:/files' cplaisier/ccafv2_seurat5
```

