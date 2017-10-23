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
data(geneList)
genes <- names(geneList)[1:100]
universe <- names(geneList)
enrichRes <- enrich.HGT(genes, universe)
enrichment_plot(enrichment=enrichRes@result, plotTitle = "Hypergemetric test")
}

