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
RankAB(beta, genelist=c(), top=10, bottom=10, scale_cutoff=1, main=NULL,filename=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
The output data frame of GroupAB or a data frame object in which the column 'Gene', ctrlname,  treatname, 'group', 'diff', and 'ENTREZID' are accessible.
}
  \item{genelist}{
An interested gene list, they will be labeled.
}
\item{top}{
     number of top genes labeled in rank figure
}
\item{bottom}{
     number of bottom genes labeled in rank figure
}
  \item{scale_cutoff}{
Cutoff scale, can be TRUE or FALSE, or a number. e.g. TRUE means the cutoff will be scaled according to library size, 1 means using one standard deviation as cutoff(default).
}
  \item{main}{
Title of rank figure
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
##Rank beta score changes and mark top and bottom 10 genes
RankAB(ddAB_essential,genelist=c(),main="Cell cycle  normalized")
}

