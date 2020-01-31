# MAGeCKFlute
This package implements methods to perform quality control (QC), normalization, batch effect removal, gene hit identification and downstream functional enrichment analysis for CRISPR screens. Before using this package, please finish the preliminary analysis using [MAGeCK](https://sourceforge.net/p/mageck/wiki/Home/) or [MAGeCK-VISPR](https://bitbucket.org/liulab/mageck-vispr/src/master/). For more detail, please read our paper [Integrative analysis pipeline for pooled CRISPR functional genetic screens](https://www.nature.com/articles/s41596-018-0113-7).


## Install package MAGeCKFlute

~~~
# Bioconductor version (stable)
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("MAGeCKFlute")

# The newest version (devel): may include new features
install.packages("devtools")
library(devtools)
install_bitbucket("liulab/MAGeCKFlute")
~~~

## Install package MAGeCKFlute using conda

~~~
conda install -c bioconda bioconductor-mageckflute
~~~


## Documentation
Details on how to use MAGeCKFlute are available on bioconductor website:
https://www.bioconductor.org/packages/3.10/bioc/vignettes/MAGeCKFlute/inst/doc/MAGeCKFlute.html
https://www.bioconductor.org/packages/3.10/bioc/vignettes/MAGeCKFlute/inst/doc/MAGeCKFlute_enrichment.html


## Version history
Please browse the [version tracking page](https://www.bioconductor.org/packages/3.10/bioc/news/MAGeCKFlute/NEWS).
	
## Contacts

* Wubing Zhang (Watson5bZhang@gmail.com)
* Binbin Wang (wangbinbintj@gmail.com)
* Wei Li (li.david.wei@gmail.com)
