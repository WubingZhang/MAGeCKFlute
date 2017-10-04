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
Square.plot(beta, ctrlname = "Control", treatname = "Treatment", main = NULL, filename = NULL, out.dir = ".", scale_cutoff = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
A data frame object in which the column "Gene","Treatment","Control","ENTREZID" are accessible.
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
  \item{scale_cutoff}{
Cutoff scale, can be TRUE or FALSE, or a number. e.g. TRUE means the cutoff will be scaled according to library size, 1 means using one standard deviation as cutoff(default).
}
}

\value{
A list contains a ggplot object and a data frame.
\item{p }{The scatter plot object}
\item{dd1 }{A data frame contains Group information of all genes}

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

dd=dd[,c("Gene","Treatment","Control","ENTREZID")]
P2=Square.plot(dd,main="Cell cycle normalized")
print(P2$p)
head(P2$dd1)
}
