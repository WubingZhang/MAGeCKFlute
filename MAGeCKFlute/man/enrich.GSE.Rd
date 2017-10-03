\name{enrich.GSE}
\alias{enrich.GSE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Do enrichment analysis using GSEA
}
\description{
Run GSEA in R easily
}
\usage{
enrich.GSE(geneList, type = "KEGG", organism="hsa",
  minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "none")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneList}{
a ranked numeric list with geneid as the names, and names of this genelist will be served as backgound.
}
  \item{type}{
gene set type, can be one of "KEGG", "BP", "MF", "CC", "DO", "NCG" and "MKEGG"
}
\item{organism}{
     "hsa" or "mmu"
}
  \item{minGSSize}{
the min size of gene set to do enrichment test
}
  \item{maxGSSize}{
the max size of gene set to do enrichment test
}
  \item{pvalueCutoff}{
the cutoff of pvalue and adjust pvalue
}
  \item{pAdjustMethod}{
p value correction method. Can be abbreviated.
}
}
\value{
A enrichResult object
}

\author{
Wubing Zhang, Feizhen Wu, Binbin Wang
}

\seealso{
 \code{\link{enrich.HGT}}
 \code{\link{enrich.GOstats}}
 \code{\link{enrich.DAVID}}
 \code{\link{enrich.HGT}}
 \code{\link{enrich.ORT}}
 \code{\link{enrichment_analysis}}
 \code{\link{enrichResult-class}}
}

\examples{
data("BRAF_mle.gene_summary")
gene_summary = BRAF_mle.gene_summary

#Read beta score from gene summary matrix
dd = ReadBeta(gene_summary, ctrlName = "D7_R1", treatName = "PLX7_R1")

##Select GroupA and GroupB genes
ddAB=GroupAB(dd)

geneList=ddAB$diff
names(geneList)=ddAB$ENTREZID

enrichRes = enrich.GSE(geneList, type = "KEGG", organism="hsa")
head(enrichRes@result)

}
