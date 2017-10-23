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
  cutoff=Cutoff_Calling(rnorm(18000),scale=TRUE)
}
