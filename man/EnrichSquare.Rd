\name{EnrichSquare}
\alias{EnrichSquare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Do enrichment analysis for selected genes in 9 Square plot
}
\description{
For each group of selected genes and combinations of them in 9 Square plot, do enrichment analysis on KEGG pathways and GO BP terms.
}
\usage{
EnrichSquare(beta, ctrlname = "Control", treatname = "Treatment",
            pvalue = 0.05, enrich_method="Hypergeometric", organism="hsa",
            adjust = "BH", filename = NULL, out.dir = ".")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
The output data frame of Square.plot or a data frame object in which the column Gene, ctrlname,  treatname, group and ENTREZID are accessible.
}
  \item{ctrlname}{
Label of Control sample
}
  \item{treatname}{
Label of Treatment sample
}
  \item{pvalue}{
Pvalue cutoff used in enrichment analysis
}
  \item{enrich_method}{
Enrichment method for KEGG pathway enrichment
}
\item{organism}{
     organism, such as ("hsa" or Human), ("mmu" or "Mouse")
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
}

\value{
A list contains all grid plots for enrichment results, and all of the return grid plots are ggplot object.
%%  If it is a LIST, use
  \item{kegg1 }{A list including a grid plot and an enrichResult object for pathway enrichmnent results of Group1 genes}
  \item{kegg2 }{A list including a grid plot and an enrichResult object for pathway enrichmnent results of Group2 genes}
  \item{kegg3 }{A list including a grid plot and an enrichResult object for pathway enrichmnent results of Group3 genes}
  \item{kegg4 }{A list including a grid plot and an enrichResult object for pathway enrichmnent results of Group4 genes}
  \item{kegg13 }{A list including a grid plot and an enrichResult object for pathway enrichmnent results of the combination of Group1 and Group3 genes}
  \item{kegg14 }{A list including a grid plot and an enrichResult object for pathway enrichmnent results of the combination of Group1 and Group4 genes}
  \item{kegg23 }{A list including a grid plot and an enrichResult object for pathway enrichmnent results of the combination of Group2 and Group3 genes}
  \item{kegg24 }{A list including a grid plot and an enrichResult object for pathway enrichmnent results of the combination of Group2 and Group4 genes}
  \item{kegg1234 }{A list including a grid plot and an enrichResult object for pathway enrichmnent results of the combination of 4 group genes}

  \item{bp1 }{A list including a grid plot and an enrichResult object including a grid plot and an enrichResult object for BP enrichmnent results of Group1 genes}
  \item{bp2 }{A list including a grid plot and an enrichResult object for BP enrichmnent results of Group2 genes}
  \item{bp3 }{A list including a grid plot and an enrichResult object for BP enrichmnent results of Group3 genes}
  \item{bp4 }{A list including a grid plot and an enrichResult object for BP enrichmnent results of Group4 genes}
  \item{bp13 }{A list including a grid plot and an enrichResult object for BP enrichmnent results of the combination of Group1 and Group3 genes}
  \item{bp14 }{A list including a grid plot and an enrichResult object for BP enrichmnent results of the combination of Group1 and Group4 genes}
  \item{bp23 }{A list including a grid plot and an enrichResult object for BP enrichmnent results of the combination of Group2 and Group3 genes}
  \item{bp24 }{A list including a grid plot and an enrichResult object for BP enrichmnent results of the combination of Group2 and Group4 genes}
  \item{bp1234 }{A list including a grid plot and an enrichResult object for BP enrichmnent results of the combination of 4 group genes}

The grid plot is a ggplot object, and Class "enrichResult" represents the result of enrichment analysis. Check more detail about enrichResult please see ?`enrichResult-class`.
}

\author{
Wubing Zhang, Feizhen Wu, Binbin Wang
}
\seealso{
\code{\link{Square.plot}}
}
\examples{
# gene_summary = BRAF_mle.gene_summary
# #Read beta score from gene summary matrix
# dd = ReadBeta(gene_summary, ctrlName = "D7_R1",
#               treatName = "PLX7_R1", organism="hsa")
# head(dd)
# dd_essential = NormalizeBeta(dd, method="cell_cycle")

# dd2=dd_essential[,c("Gene","Treatment","Control","ENTREZID")]
# P2=Square.plot(dd2,main="Cell cycle normalized")

# # Enrichment analysis for grouped genes
# E1 = EnrichSquare(P2$dd1, pvalue=0.05, adjust="BH",
#                   enrich_method="ORT", organism="hsa")
# print(E1$kegg1$gridPlot)
}

