\name{RankAB}
\alias{RankAB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rank all genes according to beta score difference between Control and Treatment samples
}
\description{
Rank all genes accounding to beta score difference between Control and Treatment samples, and label top and bottom meaningful genes. Some other interested genes can be labeled with gene symbol.
}
\usage{
RankAB(beta, genelist = c(), cutoff = Cutoff_Calling(data$diff, scale = T), main = NULL, filename = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
The output data frame of GroupAB or a data frame object in which the column 'Gene', ctrlname,  treatname, 'group', 'diff', and 'ENTREZID' are accessible.
}
  \item{genelist}{
An interested gene list, they will be labeled.
}
  \item{cutoff}{
The cutoff in use to filter GroupA and GroupB genes
}
  \item{main}{
Title of rank figure
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
##Rank beta score changes and mark top and bottom 10 genes
RankAB(ddAB_essential$result,genelist=c(),ddAB_essential$cutoff,main="Cell cycle  normalized")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
