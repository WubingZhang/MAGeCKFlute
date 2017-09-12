\name{Square.plot}
\alias{Square.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Select drug related genes using 9 Square model
}
\description{
Plot a scatter plot with Control beta score as x-axis and Treatment beta score as y-axis. Select meaningful genes by considering the detail information of Control beta score, Treatment beta score and the beta score difference  in two samples.
}
\usage{
Square.plot(beta, ctrlname = "Control", treatname = "Treatment", main = NULL, filename = NULL, out.dir = ".")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
%%     ~~Describe \code{beta} here~~
}
  \item{ctrlname}{
Label of Control sample
}
  \item{treatname}{
Label of Treatment sample
}
  \item{main}{
Title of figure
}
  \item{filename}{
The file name of output figure, when it is NULL(default), figure will not be output.
}
  \item{out.dir}{
Output directory
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list contains a ggplot object and a data frame.
\item{p }{The scatter plot object}
\item{dd1 }{A data frame contains Group information of all genes}

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
dd_essential = NormalizeBeta(dd)

dd2=dd_essential[,c("Gene","Treatment","Control","ENTREZID")]
P2=Square.plot(dd2,main="Cell cycle normalized")
print(P2$p)
head(P2$dd1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
