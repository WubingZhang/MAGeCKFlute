[![platform](https://bioconductor.org/shields/availability/3.11/MAGeCKFlute.svg)](https://bioconductor.org/shields/availability/3.11/MAGeCKFlute) [![Build Status](https://travis-ci.com/WubingZhang/MAGeCKFlute.svg)](https://travis-ci.com/WubingZhang/MAGeCKFlute) [![Update Status](https://www.bioconductor.org/shields/lastcommit/release/bioc/MAGeCKFlute.svg)](https://www.bioconductor.org/packages/3.11/bioc/html/MAGeCKFlute.html)
# MAGeCKFlute 

This package implements methods to perform quality control (QC), normalization, batch effect removal, gene hit identification and downstream functional enrichment analysis for CRISPR screens. Before using this package, please finish the preliminary analysis using [MAGeCK](https://sourceforge.net/p/mageck/wiki/Home/) or [MAGeCK-VISPR](https://bitbucket.org/liulab/mageck-vispr/src/master/). For more detail, please read our paper [Integrative analysis pipeline for pooled CRISPR functional genetic screens](https://www.nature.com/articles/s41596-018-0113-7).


## Installation
Installing the package in a fresh R environment may take a long time. It may fail because of some issues. You can check the possible issues and solutions from https://github.com/WubingZhang/MAGeCKFlute/issues/3, or post a new issue there.


### Prerequisites
To install MAGeCKFlute, you have to first install conda following the document (https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html#install-macos-silent). Afterwards you have to register the bioconda and conda-forge channels as a package source for conda.

~~~
$ conda config --add channels conda-forge
$ conda config --add channels bioconda
$ conda config --add channels r
~~~

Before start installation, you may also need to install the libxml2 which is required by the R package XML.
~~~
conda install -c anaconda libxml2
~~~

### Installation using R
~~~
> install.packages(c("devtools", "BiocManager"), repos = "https://cloud.r-project.org")
> BiocManager::install(c("pathview", "biomaRt", "msigdbr", "dendextend", "pheatmap", "sva", "ggsci", "ggrepel", "ggpubr", "knitr", "clusterProfiler"))
> devtools::install_bitbucket("liulab/MAGeCKFlute") # Development version
# Or
> BiocManager::install("MAGeCKFlute") # Released version
~~~

## Documentation
Details on how to use MAGeCKFlute are available on bioconductor website:
https://www.bioconductor.org/packages/3.11/bioc/vignettes/MAGeCKFlute/inst/doc/MAGeCKFlute.html
https://www.bioconductor.org/packages/3.11/bioc/vignettes/MAGeCKFlute/inst/doc/MAGeCKFlute_enrichment.html


## Version history
Please browse the [version tracking page](https://www.bioconductor.org/packages/3.11/bioc/news/MAGeCKFlute/NEWS).
	
## Contacts

* Wubing Zhang (Watson5bZhang@gmail.com)
* Binbin Wang (wangbinbintj@gmail.com)
* Wei Li (li.david.wei@gmail.com)
