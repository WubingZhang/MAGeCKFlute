\name{FluteRRA}
\alias{FluteRRA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functional analysis for RRA results
}
\description{
This function do a functional analysis for essential genes identified by MAGeCK-RRA. KEGG pathway and BP enrichment are done.
}
\usage{
FluteRRA(gene_summary, prefix="Test", enrich_kegg="HyperGeometric", organism="hsa", pvalueCutoff=0.05,
                     adjust="none", out.dir=paste0(prefix, "_Flute_Results"), workspace=".")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene_summary}{
    Either a path of gene summary file in RRA results or a data frame in which the column 'id', neg.fdr are accessible.
}
  \item{prefix}{
     Prefix of output folder and result summary file
}
  \item{enrich_kegg}{
     Method used for KEGG pathway enrichment analysis
}
\item{organism}{
     "hsa" or "mmu"
}
  \item{pvalueCutoff}{
     Cutoff of pvalue used to select essential genes and do enrichment analysis
}
  \item{adjust}{
     The pvalue adjust method for enrichment analysis
}
  \item{out.dir}{
     Path of output directory
}
  \item{workspace}{
     Path of workspace
}
}

\author{
Wubing Zhang, Feizhen Wu, Binbin Wang
}

\examples{
##Run the Flute RRA pipeline
FluteRRA(BRAF_rra.gene_summary, prefix="BRAF", organism="hsa")
}

