% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RankView.R
\docType{methods}
\name{RankView}
\alias{RankView}
\alias{rankview}
\title{View the rank of gene points}
\usage{
RankView(
  rankdata,
  genelist = NULL,
  top = 10,
  bottom = 10,
  cutoff = NULL,
  main = NULL,
  filename = NULL,
  width = 5,
  height = 4,
  ...
)
}
\arguments{
\item{rankdata}{Numeric vector, with gene as names.}

\item{genelist}{Character vector, specifying genes to be labeled in figure.}

\item{top}{Integer, specifying number of top genes to be labeled.}

\item{bottom}{Integer, specifying number of bottom genes to be labeled.}

\item{cutoff}{Numeric.}

\item{main}{As in 'plot'.}

\item{filename}{Figure file name to create on disk. Default filename="NULL", which means no output.}

\item{width}{As in ggsave.}

\item{height}{As in ggsave.}

\item{...}{Other available parameters in function 'ggsave'.}
}
\value{
An object created by \code{ggplot}, which can be assigned and further customized.
}
\description{
Rank all genes according to beta score deviation, and label top and bottom meaningful genes.
Some other interested genes can be labeled too.
}
\examples{
file1 = file.path(system.file("extdata", package = "MAGeCKFlute"),
"testdata/rra.gene_summary.txt")
gdata = ReadRRA(file1)
rankdata = gdata$Score
names(rankdata) = gdata$id
RankView(rankdata)

}
\author{
Wubing Zhang
}
