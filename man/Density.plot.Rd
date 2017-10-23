\name{Density.plot}
\alias{Density.plot}

\title{
Density plot
}
\description{
Density.plot() plots the density of two sample beta scores.
}
\usage{
Density.plot(beta, ctrlname = "Control", treatname = "Treatment",
              main = NULL, xlab = "Beta Score", filename = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
A data frame object in which the column Gene, ctrlname.beta and treatname.beta are accessible.
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
  \item{xlab}{
The text for x axis
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
\seealso{
  \code{\link[ggplot2]{geom_density}}
  \code{\link{Density.diff}}
  \code{\link{Violin.plot}}
  \code{\link{Square.plot}}
}
\examples{
gene_summary = BRAF_mle.gene_summary
#Read beta score from gene summary matrix
dd = ReadBeta(gene_summary, ctrlName = "D7_R1", treatName = "PLX7_R1")
##Density plot of all gene beta scores in control and treatment
Density.plot(dd_essential)

}



