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
  pvalueCutoff = 0.05, pAdjustMethod = "none")
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
     "hsa" or "mmu"
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
 \code{\link{enrichResult-class}}
}
\examples{
data("BRAF_mle.gene_summary")
gene_summary = BRAF_mle.gene_summary

#Read beta score from gene summary matrix
dd = ReadBeta(gene_summary, ctrlName = "D7_R1", treatName = "PLX7_R1")

##Select GroupA and GroupB genes
ddAB=GroupAB(dd)

##Intersted gene list
idx1=ddAB$group=="up"
genes=ddAB$ENTREZID[idx1]
universe=ddAB$ENTREZID

enrichRes = enrich.GOstats(gene=genes, universe=universe, type="BP")
head(enrichRes@result)

}
