\name{getOrg}
\alias{getOrg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get annotation data for specific organism
}
\description{
Mapping the organism to annotation package and data for geneid convertor
}
\usage{
getOrg(organism, onlyLib = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{organism}{
organism, such as ("hsa" or Human), ("mmu" or "Mouse")
}
  \item{onlyLib}{
a boolean value indicates if only return the annotation package name
}
}

\value{
\item{pkg }{the annotation package name}
\item{Symbol_Entrez }{mapping data of gene entrezid and symbol}

}

\author{
Wubing Zhang, Feizhen Wu, Binbin Wang
}

\examples{
mapping = getOrg("human", onlyLib = FALSE)$Symbol_Entrez
head(mapping)
}
