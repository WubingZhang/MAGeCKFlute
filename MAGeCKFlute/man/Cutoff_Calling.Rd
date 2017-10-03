\name{Cutoff_Calling}
\alias{Cutoff_Calling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Call cutoff value
}
\description{
Calculate beta score cutoff
}
\usage{
Cutoff_Calling(d, scale = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
A numeric vector of beta scores
}
  \item{scale}{
Can be TRUE or FALSE indicating if the cutoff will be scaled according library size. Also can be a number, indicate how many standard deviation will be used as cutoff.
}
}

\value{
A numeric value as the cutoff
}

\author{
Wubing Zhang, Feizhen Wu, Binbin Wang
}

\examples{
  data("BRAF_mle.gene_summary")
  gene_summary = BRAF_mle.gene_summary
  #Read beta score from gene summary matrix
  dd = ReadBeta(gene_summary, ctrlName = "D7_R1", treatName = "PLX7_R1")

  dd$diff=dd$Treatment - dd$Control
  cutoff=Cutoff_Calling(dd$diff,scale=TRUE)
}
