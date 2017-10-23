\name{enrichment_plot_GSE}
\alias{enrichment_plot_GSE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot results of GSEA analysis
}
\description{
Grid plot for GSEA results
}
\usage{
enrichment_plot_GSE(enrichment, plotTitle = NULL, gridColour = "blue",
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
Wubing Zhang, Feizhen Wu, Binbin Wang
}
\seealso{
  \code{\link{enrichment_plot}}
}

\examples{
data(geneList)
# enrichRes = enrich.GSE(geneList, type = "KEGG", organism="hsa")
# enrichment_plot_GSE(enrichment=enrichRes@result, plotTitle = "GSEA Analysis")
}
