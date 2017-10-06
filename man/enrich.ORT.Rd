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
          pvalueCutoff = 0.05, qvalueCutoff = 0.2, pAdjustMethod = "none",
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
     "hsa" or "mmu"
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
# data("BRAF_mle.gene_summary")
# gene_summary = BRAF_mle.gene_summary

# #Read beta score from gene summary matrix
# dd = ReadBeta(gene_summary, ctrlName = "D7_R1", treatName = "PLX7_R1")

# ##Select GroupA and GroupB genes
# ddAB=GroupAB(dd)

# ##Intersted gene list
# idx1=ddAB$group=="up"
# genes=ddAB$ENTREZID[idx1]
# universe=ddAB$ENTREZID

# enrichRes = enrich.ORT(gene=genes, universe=universe, type="BP",
# organism="hsa")
# head(enrichRes@result)
}

