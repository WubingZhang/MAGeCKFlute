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
FluteRRA(gene_summary, ctrlName="Control", treatName="Treatment", pvalueCutoff=0.05,
                     enrich_kegg="HyperGeometric",adjust="none")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene_summary}{
    Either a path of gene summary file in RRA results or a data frame in which the column 'id', {treatName}_vs_{ctrlName}.neg.fdr are accessible.
}
  \item{ctrlName}{
     Label of Control sample
}
  \item{treatName}{
     Label of Treatment sample
}
  \item{enrich_kegg}{
     Method used for KEGG pathway enrichment analysis
}
  \item{adjust}{
     The pvalue adjust method for enrichment analysis
}
  \item{pvalueCutoff}{
     Cutoff of pvalue used to select essential genes and do enrichment analysis
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Wubing Zhang
}

\examples{
##Load RRA gene summary result
data("BRAF_rra.gene_summary")
##Run the Flute RRA pipeline
FluteRRA(BRAF_rra.gene_summary, prefix="BRAF")
}

