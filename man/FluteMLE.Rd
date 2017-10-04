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
FluteMLE(gene_summary, ctrlname="Control", treatname="Treatment", prefix = "Test",scale_cutoff=1,
				adjust="none", pvalueCutoff=0.05, enrich_kegg="ORT", organism="hsa", gsea=F, loess=F, view_allpath=F, workspace=".")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene_summary}{
     path of gene summary file, or a data frame in which Gene, ctrlname.beta and treatname.beta are accessible.
}
\item{ctrlname}{
     the name of control sample in gene summary
}
\item{treatname}{
     the name of treatment sample in gene summary
}
\item{prefix}{
     prefix of output folder and result file
}
\item{scale_cutoff}{
 cutoff scale, can be TRUE or FALSE, or a number. e.g. TRUE means the cutoff will be scaled according to library size, 1 means using one standard deviation as cutoff(default).
}
\item{adjust}{
     correction method. Can be "bonferroni", "holm", "hochberg", "hommel", "BH" or its alias "fdr", "BY", or "none".
}
\item{pvalueCutoff}{
     p value cutoff used in enrichment analysis
}
\item{enrich_kegg}{
     enrichment method for KEGG, One of "ORT", "GSEA", "DAVID", "GOstats", and "HyperGeometric", or index from 1 to 5
}
\item{organism}{
     "hsa" or "mmu"
}
\item{gsea}{
     boolean, whether do GSEA analysis
}
\item{loess}{
     boolean, whether include loess analysis in pipeline
}
\item{view_allpath}{
     boolean, whether output all pathway figures
}

\item{workspace}{
     workspace path
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

\author{
  Wubing Zhang, Feizhen Wu, Binbin Wang
}

\examples{
## Run the Flute MLE pipeline
# FluteMLE(BRAF_mle.gene_summary, ctrlname=c("D7_R1","D7_R2"), treatname=c("PLX7_R1","PLX7_R2"),
#          prefix="BRAF_D7", pvalueCutoff=0.05, scale_cutoff=1.1, organism="hsa")
}
