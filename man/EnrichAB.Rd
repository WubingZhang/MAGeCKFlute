\name{EnrichAB}
\alias{EnrichAB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Enrichment analysis for GroupA and GroupB genes
}
\description{
The function GroupAB() returns positive selection and negative selection genes, termed as GroupA and GroupB genes. Using EnrichAB() to do enrichment analysis for GroupA and GroupB genes.
}
\usage{
EnrichAB(data, pvalue = 0.05, enrich_method = "Hypergeometric", adjust = "none", filename = NULL, out.dir = ".", gsea = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The output data frame of GroupAB or a data frame object in which the column 'Gene', {ctrlname},  {treatname}, 'group' and 'ENTREZID' are accessible.
}
  \item{pvalue}{
Pvalue cutoff used in enrichment analysis
}
  \item{enrich_method}{
Enrichment method for KEGG pathway enrichment
}
  \item{adjust}{
Pvalue adjust method, one of "none", "BH" and "bonferroni"
}
  \item{filename}{
Suffix of output file name, if it is NULL(default), figures will not be output.
}
  \item{out.dir}{
Output directory
}
  \item{gsea}{
A boolean value indicates if include GSEA in enrichment analysis
}
}

\value{
A list contains all grid plots for enrichment results, and all of the return grid plots are ggplot object.
%%  If it is a LIST, use
  \item{keggA }{A list including a grid plot and an enrichResult object for pathway enrichmnent results of GroupA genes}
  \item{keggB }{A list including a grid plot and an enrichResult object for pathway enrichmnent results of GroupB genes}

  \item{bpA }{A list including a grid plot and an enrichResult object including a grid plot and an enrichResult object for BP enrichmnent results of Group1 genes}
  \item{bpB }{A list including a grid plot and an enrichResult object for BP enrichmnent results of Group2 genes}

The grid plot is a ggplot object, and Class "enrichResult" represents the result of enrichment analysis. Check more detail about enrichResult please see ?`enrichResult-class`.

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Wubing Zhang
}


\seealso{
\code{\link{GroupAB}}
}
\examples{
data("BRAF_mle.gene_summary")
gene_summary = BRAF_mle.gene_summary
#Read beta score from gene summary matrix
dd = ReadBeta(gene_summary, ctrlName = "D7_R1", treatName = "PLX7_R1")
head(dd)
dd_essential = NormalizeBeta(dd, method="cell_cycle")

##Rank beta score changes between treatment and control samples, and select GroupA and GroupB genes
ddAB_essential=GroupAB(dd_essential)
##Scatter plot of GroupA and GroupB genes
ScatterAB(ddAB_essential$result,intercept=ddAB_essential$cutoff,slope=1,main="Cell cycle normalized")

##BP and KEGG enrichment analysis
enrich_result = EnrichAB(ddAB_essential$result, pvalue=0.05, gsea=F)
print(enrich_result$keggA$gridPlot)
print(enrich_result$bpA$gridPlot)
print(enrich_result$gseA$gridPlot)

}
