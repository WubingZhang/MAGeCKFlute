\name{NormalizeBeta}
\alias{NormalizeBeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normalize beta score in gene summary
}
\description{
Normalize beta scores using essential normalization or loess normalization methods.
}
\usage{
NormalizeBeta(beta, method="cell_cycle", posControl=NULL, minus=0.6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
a data frame object in which the column Gene, Control, Treatment and ENTREZID are accessible.
}
  \item{method}{
one of essential or loess
}
\item{posControl}{
    a file or character vector of essential gene list used for cell cycle normalization
}
  \item{minus}{
a useful parameter for essential normalization, it should be between 0 and 1.
}
}

\value{
a data frame has four columns, Gene, Control, Treatment, and ENTREZID.
}

\author{
Wubing Zhang, Feizhen Wu, Binbin Wang
}

\examples{
gene_summary = BRAF_mle.gene_summary

#Read beta score from gene summary matrix
dd = ReadBeta(gene_summary, ctrlName = "D7_R1", treatName = "PLX7_R1")
head(dd)

#Essential normalization
dd_essential = NormalizeBeta(dd, method="cell_cycle")
head(dd_essential)

#Optional loess normalization
dd_loess = NormalizeBeta(dd, method="loess")
head(dd_loess)
}

