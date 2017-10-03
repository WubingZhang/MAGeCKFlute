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
ScatterAB(beta, ctrlname = "Control", treatname = "Treatment", scale_cutoff = 1, main = NULL, filename = NULL)
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
  \item{scale_cutoff}{
Cutoff scale, can be TRUE or FALSE, or a number. e.g. TRUE means the cutoff will be scaled according to library size, 1 means using one standard deviation as cutoff(default).
}
  \item{main}{
Title of figure
}
  \item{filename}{
The file name of output figure, when it is NULL(default), figure will not be output.
}
}

\value{
A ggplot figure object
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

##Rank beta score changes between treatment and control samples, and select GroupA and GroupB genes
ddAB_essential=GroupAB(dd_essential)
##Scatter plot of GroupA and GroupB genes
ScatterAB(ddAB_essential$result,main="Cell cycle normalized")
}

