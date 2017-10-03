\name{FluteMLE}
\alias{FluteMLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  A pipeline of functional analysis of CRISPR screening data
}
\description{
This function do a functional analysis based on  gene summary file. In this pipeline, negative control normalized, cell cycle normalized and loess normalized(optional) datasets are analyzed simultaneously. Positive and negative selection genes are filtered as GroupA and GroupB, and drug related meaningful genes are selected. For each group of selected genes, enrichment analysis on KEGG pathways and GO BP terms is done.
}
\usage{
FluteMLE(gene_summary, ctrlname = "Control", treatname = "Treatment",
      prefix = "test", out.dir = paste0(prefix, "_Results"),
      enrich_kegg = "HyperGeometric", adjust = "none", pvalueCutoff = 0.05,
      loess = F, gsea = F, workspace = ".")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene_summary}{
     Either a path of gene summary file or a data frame in which the column 'Gene', {ctrlname}.beta and  {treatname}.beta are accessible.
}
\item{ctrlname}{
     Label of Control sample
}
\item{treatname}{
     Label of Treatment sample
}
\item{prefix}{
     Prefix of output folder and result summary file
}
\item{enrich_kegg}{
     Method used for KEGG pathway enrichment analysis
}
\item{loess}{
     A boolean value indicates whether includes pipeline analysis for loess normalized data
}
\item{adjust}{
     The pvalue adjust method for enrichment analysis
}
\item{pvalueCutoff}{
     Cutoff of pvalue used in enrichment analysis
}
\item{gsea}{
     A boolean value indicates whether includes GSEA in pipeline analysis
}
\item{out.dir}{
     Path of output directory
}
\item{workspace}{
     Path of workspace
}
}

\value{
  All of the pipeline results is output into the {out.dir}/{prefix}_Results, which includes a pdf file and many folders. The pdf file '{prefix}_Pipeline_results.pdf' is the summary of pipeline results. For each section in this pipeline, figures and useful data are outputed to corresponding subfolders.

  Distribution_of_BetaScores: Density plot and violin plot of beta scores.

  MAplot: Maplot for each normalized data.

  Linear_Fitting_of_BetaScores: Linear fitting of beta scores indicates the difference of cell cycle time between Control and Treatment samples.

  Scatter_Treat_Ctrl: Positive selection and negative selection

  Enrichment_Treat-Ctrl: Enrichment analysis for positive and negative selection genes

  Pathview_Treat_Ctrl: Pathway view for top enriched pathways

  Scatter_9Square: Using 9 Square to select drug related genes

  Enrichment_9Square: Enrichment analysis for selected genes

  Pathview_9Square: Pathway view for top enriched pathways
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Wubing Zhang
}


\examples{
##Load MLE gene summary result
data("BRAF_mle.gene_summary")
##Run the Flute MLE pipeline
FluteMLE(BRAF_mle.gene_summary, ctrlname=c("D7_R1","D7_R2"), treatname=c("PLX7_R1","PLX7_R2"),
         prefix="BRAF_D7", pvalueCutoff=0.05)
}
