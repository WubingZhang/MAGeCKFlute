\name{enrich.HGT}
\alias{enrich.HGT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Do enrichment anlaysis using HyperGeometric test
}
\description{
Just simply use hyperGeometric test to do enrichment analysis
}
\usage{
enrich.HGT(gene, universe, type = "KEGG", organism = "hsa",
          pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 2,
          maxGSSize = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene}{
interested gene list
}
  \item{universe}{
background gene list
}
  \item{type}{
category of gene sets, e.g. KEGG(default)
}
  \item{organism}{
organism, such as ("hsa" or Human), ("mmu" or "Mouse")
}
  \item{pvalueCutoff}{
p value cutoff for enrichment analysis
}
  \item{pAdjustMethod}{
p value correction method. Can be abbreviated.
}
  \item{minGSSize}{
Min gene set size
}
  \item{maxGSSize}{
Max gene set size
}
}

\value{
A enrichResult object
}

\author{
Wubing Zhang, Feizhen Wu, Binbin Wang
}

\seealso{
 \code{\link{enrich.GOstats}}
 \code{\link{enrich.DAVID}}
 \code{\link{enrich.GSE}}
 \code{\link{enrich.ORT}}
 \code{\link{enrichment_analysis}}
 \code{\link[DOSE]{enrichResult-class}}
}

\examples{
data(geneList)
genes <- names(geneList)[1:100]
universe <- names(geneList)
enrichRes <- enrich.HGT(genes, universe)
head(enrichRes@result)
}
