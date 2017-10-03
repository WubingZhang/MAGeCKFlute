\name{CellCycleFit}
\alias{CellCycleFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimate the cell cycle time difference between Control and Treatment sample by linear fitting of beta scores
}
\description{
  Assign control beta scores as x coordinate, and control and treatment beta scores as y coodinate. Fit two linear lines for two different y coodinate, one of linear lines is diagonal, and the other indicates the cell cycle difference. If the slope of second line bigger than 1, then cell cycle time of control sample is less than that of treatment sample, otherwise, more than treatment sample.
}
\usage{
CellCycleFit(beta, ctrlname = "Control", treatname = "Treatment",
            main = NULL, ylab = "Beta score", filename = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
    A data frame object in which the column Gene, ctrlname.beta and  treatname.beta are accessible.
}
  \item{ctrlname}{
    Label of control sample
}
  \item{treatname}{
     Label of treatment sample
}
  \item{main}{
    The title of figure
}
  \item{ylab}{
    The text for y axis
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
##Fitting all gene beta scores
CellCycleFit(dd_essential,ylab="Beta Score",main="Cell cycle  normalized")

}
