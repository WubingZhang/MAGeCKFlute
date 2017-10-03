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
enrichment_analysis(geneList, genes, method = 1,
                  type = c("KEGG", "BP", "MF", "CC", "DO", "NCG", "MKEGG"),
                  pvalueCutoff = 0.05, qvalueCutoff = 1, pAdjustMethod = "BH",
                  minGSSize = 2, maxGSSize = 500, plotTitle = NULL, gridColour = "blue")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneList}{
A ranked numeric list with geneid as the names, and names of this genelist will be served as backgound.
}
  \item{genes}{
ENTREZ gene id of meanigful gene list
}
  \item{method}{
One of "ORT", "GSEA", "DAVID", "GOstats", and "HyperGeometric", or number index number from 1 to 5
}
  \item{type}{
One of "KEGG", "BP", "MF", "CC", "DO", "NCG", and "MKEGG"
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
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list including a grid plot and an enrichResult object
 \item{gridPlot }{A ggplot object, grid plot of top enriched terms}
 \item{enrichRes }{A enrichResult class object, representing the enrichment analysis result. Check more detail about enrichResult class please see ?`enrichResult-class`.}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Wubing Zhang
}

\examples{
data("BRAF_mle.gene_summary")
gene_summary = BRAF_mle.gene_summary
#Read beta score from gene summary matrix
dd = ReadBeta(gene_summary, ctrlName = "D7_R1", treatName = "PLX7_R1")

##Select GroupA and GroupB genes
ddAB=GroupAB(dd)$result

##Intersted gene list
idx1=ddAB$group=="up"
#Meaningful geneset
genes=ddAB$ENTREZID[idx1]
#Background gene list, whose score is useful only if GSEA analysis is done.
geneList=ddAB$diff
names(geneList)=ddAB$ENTREZID   #background

#====GO_KEGG_enrichment=====
keggA = enrichment_analysis(geneList = geneList, genes=genes, method = "Hypergeometric",type = "KEGG", pvalueCutoff = 0.05, plotTitle="KEGG: GroupA",gridColour="#e41a1c")
bpA=enrichment_analysis(geneList = geneList, genes=genes, method = "ORT", type = "BP", pvalueCutoff = 0.05, plotTitle="BP: GroupA",gridColour="#e41a1c")

head(keggA$enrichRes@result)
print(keggA$gridPlot)
print(bpA$gridPlot)

}

