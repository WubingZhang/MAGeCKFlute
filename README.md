# MAGeCKFlute
Integrative analysis pipeline for pooled CRISPR functional genetic screens


## Install package MAGeCKFlute

~~~
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("MAGeCKFlute")

#or
install.packages("devtools")
library(devtools)
install_bitbucket("liulab/MAGeCKFlute")
~~~

## Quick start
Please look at the [tutorial page](https://www.bioconductor.org/packages/3.8/bioc/vignettes/MAGeCKFlute/inst/doc/MAGeCKFlute.html).

## Version history
	

	CHANGES IN VERSION 1.2.4
	-------------------------
	  o Fix the bug about essential gene issue for different organism.
    o Update EnrichedGeneView (Rank pathways based on adjust.pvalue).
    o Merge EnrichedView and EnrichedGSEView, which allow custom x-axis.
    o Replace EnrichGSEView in enrichment_analysis by new function.
    o Replace the old SquareView function.
    o Add a new function 'EnrichedFilter' to simplify enriched function terms by computing the Jaccard similarity coefficient between enriched terms, and removing big gene sets.
	
	CHANGES IN VERSION 1.2.3
	-------------------------

    o Prioritize FluteMLE
    o Remove warnings from vignettes


	CHANGES IN VERSION 1.2.2
	-------------------------

    o Revise ReadRRA function.
    o Add ReadsgRRA function to read sgrna_summary in MAGeCK RRA results.
    o Add sgRankView function to visualize the rank of sgRNA targeting top selected genes.
    o Remove bugs in FluteMLE and FluteRRA.
    o Add additional visualization functions in FluteRRA.
    o Add EnrichedGeneView to visualize the core enriched genes in enriched pathways.
    o Integrate WikiPathways, PID terms and EHMN pathways into enrichment functions.
    o Export all the hidden functions.


	CHANGES IN VERSION 1.2.1
	-------------------------

    o Add VolcanoView to show positive and negative selected genes.
    o Improve EnrichedView and EnrichedGSEView function.
    o Add QC figures in vignettes.
    o Merge MsigDB, KEGG, and GO genesets together and integrate them into package.


	CHANGES IN VERSION 1.1.6
	-------------------------

    o Add functions to plot figures in NatureProtocol manuscript.
    o Shorten the check time.
    o Remove null figures in pathview part.


	CHANGES IN VERSION 1.1.9
	-------------------------

    o Add parameter 'pathway_limit' / 'limit' / 'gmtpath' in FluteRRA, FluteMLE, and all enrichment functions which enable users to customize gene sets for enrichment analysis.
    o Remove DAVID and GOstats, which are not recommended.
    o Add enrichment score in enrichment results from all algorithms.


	CHANGES IN VERSION 1.1.1
	-------------------------

    o Speed up the enrichment analysis.
    o Release memory in time.


	CHANGES IN VERSION 0.99.19
	-------------------------

    o Beutify figures.
    o Remove some unnecessary dependencies.


	CHANGES IN VERSION 0.99.18
	-------------------------

    o Add HeatmapView and BatchRemove.
    o Change view distribution functions which show samples separately.
    o Change function ReadBeta to be more friendly.
    o Label top ten essential genes in SquareView.
    o Change some default parameter values to be better.

	CHANGES IN VERSION 0.99.10
	-------------------------

    o Change all plot function names ended with 'View'.
    o Decrease exported functions
    o Decrease package denpendcies
    o Revise all documents


	CHANGES IN VERSION 0.99.1
	-------------------------

    o Remove some new errors, such as error trigered by no GroupA genes.
    o Allow users to input their own essential genes to do the cell cycle normalization
    o Allow users to define the number of genes labeled in rank figure, default label top 10 and bottom 10 genes
    o Add annotation of other organisms


	CHANGES IN VERSION 0.99.0
	-------------------------

    o FluteMLE and FluteRRA are two main functions in MAGeCKFlute package. FluteMLE run pipeline from gene beta scores caculated by MAGeCK MLE, while FluteRRA run pipeline based on MAGeCK RRA results.



## Contacts

* Wubing Zhang (Watson5bZhang@gmail.com)
* Binbin Wang (wangbinbintj@gmail.com)
