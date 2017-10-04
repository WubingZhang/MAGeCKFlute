\name{enrich.DAVID}
\alias{enrich.DAVID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Do enrichment analysis using DAVID
}
\description{
An update version of DAVIDWebService to do enrichment analysis
}
\usage{
enrich.DAVID(gene, universe = NULL, david.user, idType="ENTREZ_GENE_ID", minGSSize = 2, maxGSSize = 500,
annotation = "GOTERM_BP_FAT", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene}{
interested gene list
}
  \item{universe}{
background gene list
}
  \item{david.user}{
DAVID user account
}
  \item{idType}{
gene id type for enrichment analysis in DAVID
}
  \item{minGSSize}{
min gene set limited
}
  \item{maxGSSize}{
max gene set limited
}
  \item{annotation}{
category of gene sets, e.g. GOTERM_BP_FAT(default)
}
  \item{pvalueCutoff}{
p value cutoff for enrichment analysis
}
  \item{pAdjustMethod}{
p value correction method. Can be abbreviated.
}
  \item{qvalueCutoff}{
q value cutoff for enrichment analysis
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
 \code{\link{enrich.GSE}}
 \code{\link{enrich.ORT}}
 \code{\link{enrichment_analysis}}
 \code{\link[DOSE]{enrichResult-class}}
}
\examples{

# gene_summary = BRAF_mle.gene_summary
# #Read beta score from gene summary matrix
# dd = ReadBeta(gene_summary, ctrlName = "D7_R1", treatName = "PLX7_R1")

# ##Select GroupA and GroupB genes
# ddAB=GroupAB(dd)

# ##Intersted gene list
# idx1=ddAB$group=="up"
# genes=ddAB$ENTREZID[idx1]
# universe=ddAB$ENTREZID

# enrichRes = enrich.DAVID(gene=genes, universe=universe, david.user="example@example.com")
# head(enrichRes@result)

}

