% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IncorporateDepmap.R
\docType{methods}
\name{IncorporateDepmap}
\alias{IncorporateDepmap}
\title{Incorporate Depmap screen into analysis}
\usage{
IncorporateDepmap(
  dd,
  symbol = "id",
  cell_lines = NA,
  lineages = "All",
  na.rm = FALSE
)
}
\arguments{
\item{dd}{A data frame.}

\item{symbol}{A character, specifying the column name of gene symbols in the data frame.}

\item{cell_lines}{A character vector, specifying the cell lines in Depmap to be considered.}

\item{lineages}{A character vector, specifying the lineages in Depmap to be considered.}

\item{na.rm}{Boolean, indicating whether removing NAs from the results.}
}
\value{
A data frame with Depmap column attached.
}
\description{
Incorporate Depmap screen into analysis
}
\examples{
file1 = file.path(system.file("extdata", package = "MAGeCKFlute"),
"testdata/rra.gene_summary.txt")
gdata = ReadRRA(file1)
head(gdata)
\dontrun{
  gdata = IncorporateDepmap(gdata)
  head(gdata)
}
}
\author{
Wubing Zhang
}
