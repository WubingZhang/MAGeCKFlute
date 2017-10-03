# MAGeCKFlute
A pipeline of functional analysis for CRISPR screen data

##Abstract
The purpose to do CRISPR screen is to filter essential genes, pathways and to explain the background mechanisms. We developed MAGeCKFlute to perform integrated analysis of CRISPR/Cas9 screens with/without drug treatments. The MAGeCKFlute provides several strategies to remove potential biases within read counts and beta scores. The downstream analysis for CBS and TBS with the package includes identifying essential, non- essential, and drug-associated genes, and performing biological functional analysis for these genes. The package also visualizes genes in the context of pathways to better help users explore the screening data. Collectively, MAGeCKFlute enables accurate identification of essential, non-essential, drug-targeted genes, as well as their related biological functions.

##Install package MAGeCKFlute

~~~
source("http://bioconductor.org/biocLite.R")
biocLite("MAGeCKFlute")
~~~

##Quick start

~~~
library(MAGeCKFlute)##Run pipeline from MAGeCK MLE resultsdata("BRAF_mle.gene_summary")FluteMLE(BRAF_mle.gene_summary, ctrlname=c("D7_R1", "D7_R2"), treatname=c("PLX7_R1","PLX7_R2"), 
		prefix="BRAF", organism =”hsa”)

##Run pipeline from MAGeCK RRA resultsdata("BRAF_rra.gene_summary")FluteRRA(BRAF_rra.gene_summary, prefix="BRAF", organism=”hsa”)
~~~

##Contacts

* Wubing Zhang (Watson5bZhang@gmail.com)
* Binbin Wang (wangbinbintj@gmail.com)
* Feizhen Wu
