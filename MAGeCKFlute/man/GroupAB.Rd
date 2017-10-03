\name{GroupAB}
\alias{GroupAB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Negative selection and positive selection
}
\description{
Compare Control beta scores and Treatment beta scores, and then get positive selection and negative selection genes.
}
\usage{
GroupAB(beta, ctrlname = "Control", treatname = "Treatment", scale_cutoff = 1, filename = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
A data frame object in which the column 'Gene', ctrlname.beta and  treatname.beta are accessible.
}
  \item{ctrlname}{
Label of Control sample
}
  \item{treatname}{
Label of Treatment sample
}
  \item{scale_cutoff}{
Cutoff scale, can be TRUE or FALSE, or a number. e.g. TRUE means the cutoff will be scaled according to library size, 1 means using one standard deviation as cutoff(default).
}
  \item{filename}{
The file name of output figure, when it is NULL(default), figure will not be output.
}
}

\value{
A list contains dd1 object and cutoff object
  \item{dd1 }{A data frame which has column of Gene, treatname, ctrlname, ENTREZID, diff, and group}
  \item{cutoff }{The cutoff of beta score to filter negative selection and positive selection genes}

}

\author{
Wubing Zhang, Feizhen Wu, Binbin Wang
}

\examples{
data("BRAF_mle.gene_summary")
gene_summary = BRAF_mle.gene_summary
ctrlName = c("D7_R1", "D7_R2")
treatName = c("PLX7_R1", "PLX7_R2")

#Read beta score from gene summary matrix
dd = ReadBeta(gene_summary, ctrlName = ctrlName, treatName = treatName)

#Essential normalization
dd_essential = NormalizeBeta(dd, method="cell_cycle")

ddAB_essential=GroupAB(dd_essential)

##View the head of data
head(ddAB_essential)
}

