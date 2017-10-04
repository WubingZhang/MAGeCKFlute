\name{enrichment_analysis}
\alias{enrichment_analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Do enrichment analyis for gene list
}
\description{
Do enrichment analysis on KEGG pathways, BP, CC, MF, DO, NCG, or MKEGG, and various methods are available, including "ORT", "GSEA", "DAVID", "GOstats", and "HyperGeometric".
}
\usage{
enrichment_analysis(geneList, universe, method = 1, type = "KEGG", organism="hsa",
      pvalueCutoff = 0.05, qvalueCutoff = 1, pAdjustMethod = "none",
      minGSSize = 2, maxGSSize = 500, plotTitle = NULL, gridColour = "blue")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneList}{
a vector of interested genes symbols.
}
  \item{universe}{
a vector of background gene symbols
}
  \item{method}{
One of "ORT", "GSEA", "DAVID", "GOstats", and "HyperGeometric", or index from 1 to 5
}
  \item{type}{
One of "KEGG", "BP", "MF", "CC", "DO", "NCG", and "MKEGG"
}
\item{organism}{
     "hsa" or "mmu"
}
  \item{pvalueCutoff}{
The cutoff of pvalue and adjust pvalue
}
  \item{qvalueCutoff}{
The cutoff of qvalue
}
  \item{pAdjustMethod}{
One of "none", "BH" and "bonferroni"
}
  \item{minGSSize}{
The min size of gene set to do enrichment test
}
  \item{maxGSSize}{
The max size of gene set to do enrichment test
}
  \item{plotTitle}{
The title of grid plot
}
  \item{gridColour}{
The grid color of grid plot
}
}

\value{
A list including a grid plot and an enrichResult object
 \item{gridPlot }{A ggplot object, grid plot of top enriched terms}
 \item{enrichRes }{A enrichResult class object, representing the enrichment analysis result. Check more detail about enrichResult class please see ?`enrichResult-class`.}

}

\author{
Wubing Zhang, Feizhen Wu, Binbin Wang
}

\examples{

gene_summary = BRAF_mle.gene_summary
#Read beta score from gene summary matrix
dd = ReadBeta(gene_summary, ctrlName = "D7_R1", treatName = "PLX7_R1")

##Select GroupA and GroupB genes
ddAB=GroupAB(dd)

##Intersted gene list
idx1=ddAB$group=="up"
#Meaningful geneset
genes=ddAB$ENTREZID[idx1]
universe=ddAB$ENTREZID

#====GO_KEGG_enrichment=====
keggA = enrichment_analysis(geneList = genes, universe=universe, method = "Hypergeometric",type = "KEGG", organism = "hsa", plotTitle="KEGG: GroupA",gridColour="#e41a1c")

head(keggA$enrichRes@result)
print(keggA$gridPlot)

}

