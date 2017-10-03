\name{Violin.plot}
\alias{Violin.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Violin plot
}
\description{
Plots the violin of two sample beta scores.
}
\usage{
Violin.plot(beta, ctrlname = "Control", treatname = "Treatment", main = NULL, ylab = "Beta Score", filename = NULL)
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
  \item{ylab}{
The text for y axis
}
  \item{filename}{
The file name of output figure, when it is NULL(default), figure will not be output.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A ggplot figure object
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
dd_essential = NormalizeBeta(dd, method="cell_cycle")
##Violin plot of all gene beta scores in control and treatment
Violin.plot(dd_essential, main="Cell cycle normalized")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
