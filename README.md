# MAGeCKFlute 

This package implements methods to perform quality control (QC), normalization, batch effect removal, gene hit identification and downstream functional enrichment analysis for CRISPR screens. Before using this package, please finish the preliminary analysis using [MAGeCK](https://sourceforge.net/p/mageck/wiki/Home/) or [MAGeCK-VISPR](https://bitbucket.org/liulab/mageck-vispr/src/master/). For more detail, please read our paper [Integrative analysis pipeline for pooled CRISPR functional genetic screens](https://www.nature.com/articles/s41596-018-0113-7).


## Installation
Installing the package in a fresh R environment may take a long time. For any questions or comments, please post it to the [MAGeCK Google group](https://groups.google.com/d/forum/mageck).


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
> BiocManager::install(c("pathview", "biomaRt", "msigdbr", "dendextend", "pheatmap", "sva", "ggrepel", "knitr", "clusterProfiler", "depmap"))
> BiocManager::install("MAGeCKFlute") # Released version
# Or
> devtools::install_github("liulab-dfci/MAGeCKFlute")
~~~

## Documentation
Details on how to use MAGeCKFlute are available on bioconductor website: https://www.bioconductor.org/packages/release/bioc/html/MAGeCKFlute.html


## Contacts

* Binbin Wang (wangbinbintj@gmail.com)
* Wei Li (li.david.wei@gmail.com)
