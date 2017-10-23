\name{enrich.GOstats}
\alias{enrich.GOstats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Do enrichment analysis using GOstats
}
\description{
GOstats can be used easily.
}
\usage{
enrich.GOstats(gene, universe = NULL, type = c("KEGG", "BP", "MF", "CC"), organism="hsa",
  pvalueCutoff = 0.05, pAdjustMethod = "BH")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene}{
interested ENTREZID gene list
}
  \item{universe}{
background ENTREZID gene list
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
}

\value{
A enrichResult object
}

\author{
Wubing Zhang, Feizhen Wu, Binbin Wang
}


\seealso{
 \code{\link{enrich.HGT}}
 \code{\link{enrich.DAVID}}
 \code{\link{enrich.GSE}}
 \code{\link{enrich.HGT}}
 \code{\link{enrich.ORT}}
 \code{\link{enrichment_analysis}}
 \code{\link[DOSE]{enrichResult-class}}
}
\examples{
# data(geneList)
# genes <- names(geneList)[1:100]
# universe <- names(geneList)
# enrichRes <- enrich.GOstats(genes, universe, type="BP")
# head(enrichRes@result)
}
