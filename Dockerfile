FROM satijalab/seurat:5.0.0

MAINTAINER Plaisier Lab <plaisier@asu.edu>

RUN \
  apt-get update && apt-get install -y apt-transport-https && \
  apt-get install -y python3 python3-dev python3-pip python3-virtualenv libssl-dev libjpeg-dev libmagick++-dev libsodium-dev && \
  rm -rf /var/lib/apt/lists/*

RUN apt-get update -y
RUN apt-get install -y libharfbuzz-dev libfribidi-dev
RUN apt-get install -y build-essential cmake

RUN R -e "install.packages(c('remotes','roxygen2','Seurat','ggplot2'))"
RUN R -e "install.packages(c('plumber', 'yaml', 'base64enc', 'readr', 'writexl'))"
RUN R -e "install.packages(c('BiocManager'))"
RUN R -e "BiocManager::install('org.Hs.eg.db')"
RUN R -e "BiocManager::install(c('org.Mm.eg.db'))"
RUN R -e "remotes::install_github('plaisier-lab/ccafv2_R/ccAFv2')"
RUN apt-get install -y vim
RUN R -e "install.packages(c('ADImpute','yarrr'))"
