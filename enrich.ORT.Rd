\name{enrich.ORT}
\alias{enrich.ORT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Enrichment analysis using Over-Representation Test
}
\description{
Over-Representation Analysis in ClusterProfiler
}
\usage{
enrich.ORT(gene, universe = NULL, type = "KEGG", organism="hsa",
          pvalueCutoff = 0.05, qvalueCutoff = 0.2, pAdjustMethod = "BH",
          minGSSize = 2, maxGSSize = 500)
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
category of gene sets, one of "KEGG", "BP", "MF", "CC", "DO", "NCG",
and "MKEGG"
}
\item{organism}{
     organism, such as ("hsa" or Human), ("mmu" or "Mouse")
}
  \item{pvalueCutoff}{
p value cutoff for enrichment analysis
}
  \item{qvalueCutoff}{
q value cutoff
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
 \code{\link{enrich.HGT}}
 \code{\link{enrichment_analysis}}
 \code{\link[DOSE]{enrichResult-class}}
}


\examples{
data(geneList)
genes <- names(geneList)[1:100]
universe <- names(geneList)
# enrichRes <- enrich.ORT(genes, universe)
# head(enrichRes@result)
}

