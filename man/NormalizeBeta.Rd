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
NormalizeBeta(beta, method = "cell_cycle", minus = 0.6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
A data frame object in which the column Gene, Control, Treatment and ENTREZID are accessible.
}
  \item{method}{
One of essential or loess
}
  \item{minus}{
A useful parameter for essential normalization, it should be between 0 and 1.
}
}

\value{
A data frame has four columns, Gene, Control, Treatment, and ENTREZID.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Wubing Zhang
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data("BRAF_mle.gene_summary")
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
