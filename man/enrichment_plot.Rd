\name{enrichment_plot}
\alias{enrichment_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot enrichment analysis results
}
\description{
Grid plot for enrichment results
}
\usage{
enrichment_plot(enrichment, plotTitle = NULL, gridColour = "blue",
                termNum = 10, charLength = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{enrichment}{
a data frame of enrichment result
}
  \item{plotTitle}{
title of figure
}
  \item{gridColour}{
grid color in figure
}
  \item{termNum}{
term number to show
}
  \item{charLength}{
max character length to show
}
}

\value{
A ggplot object
}

\author{
Wubing Zhang, Feizhen Wu, Binbin Wu
}

\seealso{
  \code{\link{enrichment_plot_GSE}}
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

  enrichRes = enrich.HGT(gene=genes, universe=universe, type="KEGG")
  head(enrichRes@result)

  enrichment_plot(enrichment=enrichRes@result,
                  plotTitle = "Hypergemetric test")
}

