# MAGeCKFlute
Integrative analysis of pooled CRISPR functional genetic screens  

##Abstract
The purpose to do CRISPR screen is to filter essential genes, pathways and to explain the background mechanisms. We developed MAGeCKFlute to perform integrated analysis of CRISPR/Cas9 screens with/without drug treatments. The MAGeCKFlute provides several strategies to remove potential biases within read counts and beta scores. The downstream analysis for CBS and TBS with the package includes identifying essential, non- essential, and drug-associated genes, and performing biological functional analysis for these genes. The package also visualizes genes in the context of pathways to better help users explore the screening data. Collectively, MAGeCKFlute enables accurate identification of essential, non-essential, drug-targeted genes, as well as their related biological functions.


##Install MAGeCKFlute

~~~
source('http://www.bioconductor.org/biocLite.R')
biocLite('MAGeCKFlute’)
~~~

##Quick start

~~~
library(MAGeCKFlute)

##Run MAGeCKFlute with MAGeCK MLE gene summary result
data("BRAF_mle.gene_summary")
FluteMLE(gene_summary= BRAF_mle.gene_summary, ctrlname=c("D7_R1", "D7_R2"), treatname=c("PLX7_R1","PLX7_R2"), 
prefix="BRAF", organism =”hsa”)


##Run MAGeCKFlute with MAGeCK RRA gene summary result
data("BRAF_rra.gene_summary")
FluteRRA(BRAF_rra.gene_summary, prefix="BRAF", organism=”hsa”)
~~~


##Contacts
* Wubing Zhang (Watson5bZhang@gmail.com)
* Binbin Wang (wangbinbintj@gmail.com)
* Feizhen Wu