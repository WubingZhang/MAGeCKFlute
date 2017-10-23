\name{MA.plot}
\alias{MA.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MAplot of beta score
}
\description{
MAplot of beta scores in Control and Treatment samples
}
\usage{
MA.plot(beta, ctrlname = "Control", treatname = "Treatment",
        main = "Negative control normalized",
        subset = sample(1:length(M), min(c(10000, length(M)))),
        show.statistics = TRUE, span = 2/3, family.loess = "gaussian",
        cex = 1, cex.lab = 1.2, cex.axis = 1, cex.main = 1.2,
        plot.method = c("normal", "smoothScatter", "add"),
        add.loess = TRUE, lwd = 1, lty = 1, loess.col = "red",
        filename = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
A data frame object in which the column 'Gene', {ctrlname}.beta and  {treatname}.beta are accessible.
}
  \item{ctrlname}{
Label of Control sample
}
  \item{treatname}{
Label of Treatment sample
}
  \item{main}{
The title of figure
}
  \item{subset}{
The same as ma.plot
}
  \item{show.statistics}{
The same as ma.plot
}
  \item{span}{
The same as ma.plot
}
  \item{family.loess}{
The same as ma.plot
}
  \item{cex}{
The same as ma.plot
}
  \item{cex.lab}{
The same as ma.plot
}
  \item{cex.axis}{
The same as ma.plot
}
  \item{cex.main}{
The same as ma.plot
}
  \item{plot.method}{
The same as ma.plot
}
  \item{add.loess}{
The same as ma.plot
}
  \item{lwd}{
The same as ma.plot
}
  \item{lty}{
The same as ma.plot
}
  \item{loess.col}{
The same as ma.plot
}
  \item{filename}{
The same as ma.plot
}
}
\value{
No return value, just print a MAplot
}

\author{
Wubing Zhang, Feizhen Wu, Binbin Wang
}

\examples{
data("BRAF_mle.gene_summary")
gene_summary = BRAF_mle.gene_summary
#Read beta score from gene summary matrix
dd = ReadBeta(gene_summary, ctrlName = "D7_R1", treatName = "PLX7_R1")
head(dd)
dd_essential = NormalizeBeta(dd, method="cell_cycle")
MA.plot(dd_essential, cex=1, main="Cell cycle normalized")
}
