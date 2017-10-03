\name{TransGeneID}
\alias{TransGeneID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transfer genes from one type to another
}
\description{
Use packages named org.*.eg.db to transfer gene id. Now only genes of hsa and mmu are accessable.
}
\usage{
TransGeneID(genes, fromType="SYMBOL", toType="ENTREZID", organism="hsa")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genes}{
vector of genes
}
  \item{fromType}{
a character, "SYMBOL" or "ENTREZID"
}
  \item{toType}{
a character, "SYMBOL" or "ENTREZID"
}
  \item{organism}{
one of "hsa"(default) and "mmu"
}
}

\value{
A vector of return gene vector whose names are given gene names
}

\author{
Wubing Zhang, Feizhen Wu, Binbin Wang
}

\examples{
genes = BRAF_mle.gene_summary$Gene[1:10]
TransGeneID(genes, fromType="SYMBOL", toType="ENTREZID", organism="hsa")
}

