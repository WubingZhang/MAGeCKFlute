\name{pathPng2Pdf}
\alias{pathPng2Pdf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read local pathway figure into pdf file
}
\description{
Kegg pathway view and write into pdf file
}
\usage{
pathPng2Pdf(genelist, pathways = c(), organism = "hsa", view_allpath=F, title = "Group A",
  sub = "Negative control normalized", output = ".", path.archive = ".", kegg.native = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genelist}{
      a data frame with gene EntrezID as rownames, whose the column is corresponding to gene score in Control and Treatment samples.
}
  \item{pathways}{
      character vector, the KEGG pathway ID(s)
}
\item{organism}{
     organism, such as ("hsa" or Human), ("mmu" or "Mouse")
}
\item{view_allpath}{
     boolean, whether output all pathway figures
}
  \item{title}{
    title of figures
}
  \item{sub}{
    sub title of figures
}
  \item{output}{
    output path
}
  \item{path.archive}{
character, the directory of KEGG pathway data file (.xml) and image file (.png). Users may supply their own data files in the same format and naming convention of KEGG's (species code + pathway id, e.g. hsa04110.xml, hsa04110.png etc) in this directory. Default kegg.dir="." (current working directory).
}
  \item{kegg.native}{
logical, whether to render pathway graph as native KEGG graph (.png) or using graphviz layout engine (.pdf). Default kegg.native=TRUE.
}
}

\author{
Wubing Zhang, Feizhen Wu, Binbin Wang
}
\seealso{
  \code{\link{pathview}}
  \code{\link{KeggPathwayView}}
}
\examples{
gene_summary = BRAF_mle.gene_summary
#Read beta score from gene summary matrix
# dd = ReadBeta(gene_summary, ctrlName = "D7_R1", treatName = "PLX7_R1", organism="hsa")
# pathPng2Pdf(dd, pathways_hsa$PathwayID[1:2], title=NULL, sub=NULL, organism="hsa")

}
