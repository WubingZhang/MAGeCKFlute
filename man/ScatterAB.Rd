\name{ScatterAB}
\alias{ScatterAB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A scatter plot for GroupA and GroupB genes
}
\description{
A scatter plot with Control beta scores as x axis and Treatment beta scores as y axis. In this figure, GroupA genes and GroupB genes are selected and marked with diffrent color.
}
\usage{
ScatterAB(beta, ctrlname = "Control", treatname = "Treatment", intercept = Cutoff_Calling(beta[, treatname] - beta[, ctrlname]), slope = 1, main = NULL, filename = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
The output data frame of GroupAB or a data frame object in which the column 'Gene', {ctrlname},  {treatname}, 'group' and 'ENTREZID' are accessible.
}
  \item{ctrlname}{
Label of control sample
}
  \item{treatname}{
Label of Treatment sample
}
  \item{intercept}{
The intercept of slash cutoff
}
  \item{slope}{
The slope of slash cutoff
}
  \item{main}{
Title of figure
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

##Rank beta score changes between treatment and control samples, and select GroupA and GroupB genes
ddAB_essential=GroupAB(dd_essential)
##Scatter plot of GroupA and GroupB genes
ScatterAB(ddAB_essential$result,intercept=ddAB_essential$cutoff,slope=1,main="Cell cycle normalized")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
