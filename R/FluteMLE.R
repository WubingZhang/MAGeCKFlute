#' Downstream analysis based on MAGeCK-MLE result
#'
#' Integrative analysis pipeline using the gene summary table in MAGeCK MLE results
#'
#' @docType methods
#' @name FluteMLE
#' @rdname FluteMLE
#' @aliases flutemle
#'
#' @param gene_summary A data frame, which contains columns of 'Gene',
#' \code{ctrlname}.beta and \code{treatname}.beta.
#' @param ctrlname A character vector, specifying the names of control samples.
#' @param treatname A character vector, specifying the names of treatment samples.
#' @param keytype "Entrez" or "Symbol".
#' @param organism "hsa" or "mmu".
#'
#' @param scale_cutoff Boolean or numeric, whether scale cutoff to whole genome level,
#' or how many standard deviation will be used as cutoff.
#' @param top An integer, specifying number of top selected genes to be labeled in rank figure.
#' @param bottom An integer, specifying number of bottom selected genes to be labeled in rank figure.
#' @param interestGenes A character vector, specifying interested genes to be labeled in rank figure.
#'
#' @param limit A two-length vector, specifying the minimal and
#' maximal size of gene sets for enrichent analysis.
#' @param pvalueCutoff A numeric, specifying pvalue cutoff of enrichment analysis, default 1.
#' @param enrich_kegg One of "ORT"(Over-Representing Test), "GSEA"(Gene Set Enrichment Analysis), and "HGT"(HyperGemetric test).
#'
#' @param posControl A character vector, specifying a list of positive control gene symbols.
#' @param loess Boolean, whether include loess normalization in the pipeline.
#'
#' @param prefix A character, indicating the prefix of output file name, which can't contain special characters.
#' @param width The width of summary pdf in inches.
#' @param height The height of summary pdf in inches.
#' @param outdir Output directory on disk.
#' @param view_allpath Boolean, whether output all pathway view figures.
#'
#' @author Wubing Zhang
#'
#' @return All of the pipeline results is output into the \code{out.dir}/\code{prefix}_Results,
#' which includes a pdf file and many folders. The pdf file '{prefix}_Pipeline_results.pdf' is the
#' summary of pipeline results. For each section in this pipeline, figures and useful data are
#' outputed to corresponding subfolders.
#' \itemize{
#'   \item {Distribution_of_BetaScores}: {Density plot and violin plot of beta scores.}
#'   \item {MAplot}: {Maplot for each normalized data.}
#'   \item {Linear_Fitting_of_BetaScores}: {Linear fitting of beta scores indicates the difference of cell cycle
#'   time between Control and Treatment samples.}
#'   \item {Scatter_Treat_Ctrl}: {Positive selection and negative selection.}
#'   \item {Enrichment_Treat-Ctrl}: {Enrichment analysis for positive and negative selection genes.}
#'   \item {Pathview_Treat_Ctrl}: {Pathway view for top enriched pathways.}
#'   \item {Scatter_9Square}: {Using 9 Square to select drug related genes.}
#'   \item {Enrichment_9Square}: {Enrichment analysis for selected genes.}
#'   \item {Pathview_9Square}: {Pathway view for top enriched pathways.}
#' }
#'
#' @details MAGeCK-MLE can be used to analyze screen data from multi-conditioned experiments. MAGeCK-MLE
#' also normalizes the data across multiple samples, making them comparable to each other. The most important
#' ouput of MAGeCK MLE is `gene_summary` file, which includes the beta scores of multiple conditions and
#' the associated statistics. The `beta score`  for each gene describes how the gene is selected: a positive
#' beta score indicates a positive selection, and a negative beta score indicates a negative selection.
#'
#' The downstream analysis includes identifying essential, non-essential, and target-associated genes,
#' and performing biological functional category analysis and pathway enrichment analysis of these genes.
#' The function also visualizes genes in the context of pathways to benefit users exploring screening data.
#'
#'
#' @seealso \code{\link{FluteRRA}}
#'
#' @examples
#' data(mle.gene_summary)
#' \dontrun{
#'   # functional analysis for MAGeCK MLE results
#'   FluteMLE(mle.gene_summary, ctrlname = "dmso", treatname = "plx",
#'            prefix = "PLX", pvalueCutoff = 0.25, organism = "hsa")
#' }
#'
#' @import ggplot2 stats grDevices utils gridExtra grid
#' @export

FluteMLE <- function(gene_summary, ctrlname, treatname,
                     keytype = "Symbol", organism = "hsa", # Input dataset
                     norm_method = "none", posControls = NULL,
                     top = 10, bottom = 10, toplabels = NA,
                     scale_cutoff = 2, limit = c(2, 200),
                     pvalueCutoff=0.25, prefix = "",
                     width = 10, height = 7, outdir = "."){

	## Prepare the running environment ##
  {
    message(Sys.time(), " # Create output dir and pdf file...")
    outdir = file.path(outdir, paste0(prefix, "_Flute_Results"))
    dir.create(file.path(outdir), showWarnings = FALSE)
    output_pdf = paste0(prefix, "_Flute.mle_summary.pdf")
    pdf(file.path(outdir, output_pdf), width = width, height = height)
  }

  ## Beta Score Preparation ##
  {
    beta = ReadBeta(gene_summary)
    if(keytype == "Symbol")
      beta$EntrezID = TransGeneID(beta$Gene, "Symbol", "Entrez", organism = organism)
    else{
      beta$EntrezID = beta$Gene
      beta$Gene = TransGeneID(beta$Gene, "Entrez", "Symbol", organism = organism)
    }
    idx1 = is.na(beta$EntrezID)
    idx2 = !is.na(beta$EntrezID) & duplicated(beta$EntrezID)
    idx = idx1|idx2
    if(sum(idx1)>0) warning(sum(idx1), " genes are not eligible: ",
                           paste0(beta$Gene[idx1], collapse = ", "))
    if(sum(idx2)>0) warning(sum(idx2), " genes have duplicate entrez id: ",
                            paste0(beta$Gene[idx2], collapse = ", "))
    # beta = beta[!idx, ]

    if(!all(c(ctrlname, treatname) %in% colnames(beta)))
      stop("Sample name doesn't match !!!")

    dd = beta[, c("Gene", "EntrezID", ctrlname, treatname)]
    if(tolower(norm_method)=="posctrl")
      dd = NormalizeBeta(dd, samples = c(ctrlname, treatname),
                         method = "cell_cycle", posControl = posControl)
    if(tolower(norm_method)=="loess")
      dd = NormalizeBeta(dd, samples = c(ctrlname, treatname), method = "loess")
  }

	## Distribution of all genes ##
	{
	  outputDir1 = file.path(outdir, "Distribution_of_BetaScores")
	  dir.create(outputDir1, showWarnings = FALSE)
	  idx_distr = c(ctrlname, treatname)
	  P1 = ViolinView(dd[, idx_distr], ylab = "Beta score",
	                  filename = file.path(outputDir1, "violin_all.png"))
	  P2 = DensityView(dd[, idx_distr], xlab = "Beta score",
	                   filename = file.path(outputDir1, "density_all.png"))
	  P3 = ConsistencyView(dd, ctrlname, treatname, ylab = "Beta score",
	                       filename = file.path(outputDir1, "scatter_all.png"))
	  grid.arrange(P1, P2, ncol = 2)

	  ## Distribution of essential genes ##
	  data(Zuber_Essential)
	  if(is.null(posControl))
	    idx = toupper(dd$Gene) %in% toupper(Zuber_Essential$GeneSymbol)
	  else
	    idx = toupper(dd$Gene) %in% toupper(posControl)

	  if(sum(idx) > 6){
	    outputDir3 = file.path(outdir, "Linear_Fitting_of_BetaScores")
	    dir.create(outputDir3, showWarnings = FALSE)

	    #Negative control normalized
	    P1 = ViolinView(dd[idx, idx_distr], ylab = "Essential.B.S.",
	                    filename = file.path(outputDir1, "violin_ess.png"))
	    P2 = DensityView(dd[idx, idx_distr], xlab = "Essential.B.S.",
	                     filename = file.path(outputDir1, "density_ess.png"))
	    P3 = ConsistencyView(dd[idx, idx_distr], ctrlname, treatname,
	                         ylab = "Essential.B.S.",
	                         filename = file.path(outputDir1, "scatter_ess.png"))
	    grid.arrange(P1, P2, P3, ncol = 2)
	  }
	}

  ## Combine replicates ##
  {
    dd$Control = rowMeans(dd[, ctrlname, drop = FALSE])
    dd$Treatment = rowMeans(dd[, treatname, drop = FALSE])
    dd.diff = dd$Treatment - dd$Control
    names(dd.diff) = dd$Gene
  }

	## Drug-targeted genes ##
	{
	  outputDir4 = file.path(outdir, "Scatter_Treat_Ctrl")
	  dir.create(outputDir4, showWarnings = FALSE)

	  # Negative normalized
	  P1 = ScatterView(dd, main = "Negative control normalized", scale_cutoff = scale_cutoff,
	                filename = file.path(outputDir4,"Scatter_Treat-Ctrl.png"))
	  P2 = RankView(dd.diff, genelist = interestGenes, top = top, bottom = bottom, main = "Negative control normalized",
	              cutoff = CutoffCalling(dd.diff, scale = scale_cutoff),
	              filename = file.path(outputDir4, "Rank_Treat-Ctrl.png"))
	  grid.arrange(P1, P3, P2, P4, ncol = 2)
	}

	## Enrichment analysis of negative and positive selected genes ##
	{
	  outputDir5 = file.path(outdir, "Enrichment_Treat-Ctrl/")
	  outputDir6 = file.path(outdir, "Pathview_Treat_Ctrl/")
	  dir.create(outputDir5, showWarnings=FALSE)
	  dir.create(outputDir6, showWarnings=FALSE)

	  E1 = EnrichAB(P1$data, pvalue = pvalueCutoff, enrich_method = enrich_kegg,
	                organism = organism, limit = limit,
	                filename = "Negative_ctrl_normalized", out.dir = outputDir5)
	  # EnrichedView
	  grid.arrange(E1$keggA$gridPlot, E1$goA$gridPlot,
	               E1$keggB$gridPlot, E1$goB$gridPlot,
	               ncol = 2, top = "Negative control normalized")

	  # Pathway view for top 4 pathway
	  if(!is.null(E1$keggA$enrichRes) && nrow(E1$keggA$enrichRes@result)>0)
	    arrangePathview(dd, gsub("KEGG_", "", E1$keggA$enrichRes@result$ID), top = 4, ncol = 2,
	                    title="Group A", sub="Negative control normalized",organism=organism,
	                    view_allpath=view_allpath, output=outputDir6)
	  if(!is.null(E1$keggB$enrichRes) && nrow(E1$keggB$enrichRes@result)>0)
	    arrangePathview(dd, gsub("KEGG_", "", E1$keggB$enrichRes@result$ID), top = 4, ncol = 2,
	                    title="Group B", sub="Negative control normalized", organism=organism,
	                    view_allpath=view_allpath, output=outputDir6)
  }
  suppressWarnings(rm(E1, E2, E3, outputDir5, outputDir6))

	## Nine-squares ##
	{
	  dir.create(file.path(outdir, "Scatter_9Square"), showWarnings=FALSE)
	  outputDir6 <- file.path(outdir, "Scatter_9Square/Square9_")
	  P1 = SquareView(dd, label="Gene", main="Negative control normalized",
	                  x_cutoff = CutoffCalling(dd$Control, scale_cutoff),
	                  y_cutoff = CutoffCalling(dd$Treatment, scale_cutoff),
	                  intercept = CutoffCalling(dd$Treatment-dd$Control, scale_cutoff),
	                  filename=paste0(outputDir6,"scatter.png"))
	  grid.arrange(P1, ncol = 1)
	}

	## Nine-Square grouped gene enrichment ##
  outputDir7 = file.path(outdir, "Enrichment_9Square")
  outputDir8 = file.path(outdir, "Pathview_9Square")
  dir.create(outputDir7, showWarnings=FALSE)
  dir.create(outputDir8, showWarnings=FALSE)
	{# Negative control normalization
	  E1 = EnrichSquare(P1$data, pvalue = pvalueCutoff,
	                    enrich_method = enrich_kegg, organism=organism,
	                    filename="negative_normalized", limit = limit,
	                    out.dir=outputDir7)
    # EnrichView
	  grid.arrange(E1$kegg1$gridPlot, E1$go1$gridPlot, ncol = 1)
	  grid.arrange(E1$kegg2$gridPlot, E1$go2$gridPlot, ncol = 1)
	  grid.arrange(E1$kegg3$gridPlot, E1$go3$gridPlot, ncol = 1)
	  grid.arrange(E1$kegg4$gridPlot, E1$go4$gridPlot, ncol = 1)
	  grid.arrange(E1$kegg12$gridPlot, E1$go12$gridPlot, ncol = 1)
	  grid.arrange(E1$kegg13$gridPlot, E1$go13$gridPlot, ncol = 1)
	  grid.arrange(E1$kegg24$gridPlot, E1$go24$gridPlot, ncol = 1)
	  grid.arrange(E1$kegg34$gridPlot, E1$go34$gridPlot, ncol = 1)

	  # PathwayView
	  if(!is.null(E1$kegg1$enrichRes) && nrow(E1$kegg1$enrichRes@result)>0)
	    arrangePathview(dd, gsub("KEGG_", "", E1$kegg1$enrichRes@result$ID), ncol = 2,
	                    title = "Group 1", sub = "Negative control normalized",
	                    organism=organism, view_allpath=view_allpath, output=outputDir8)
	  if(!is.null(E1$kegg2$enrichRes) && nrow(E1$kegg2$enrichRes@result)>0)
	    arrangePathview(dd, gsub("KEGG_", "", E1$kegg2$enrichRes@result$ID), ncol = 2,
	                    title = "Group 2", sub = "Negative control normalized",
	                    organism=organism,view_allpath=view_allpath, output=outputDir8)
	  if(!is.null(E1$kegg3$enrichRes) && nrow(E1$kegg3$enrichRes@result)>0)
	    arrangePathview(dd, gsub("KEGG_", "", E1$kegg3$enrichRes@result$ID), ncol = 2,
	                    title = "Group 3", sub = "Negative control normalized",
	                    organism=organism, view_allpath=view_allpath, output=outputDir8)
	  if(!is.null(E1$kegg4$enrichRes) && nrow(E1$kegg4$enrichRes@result)>0)
	    arrangePathview(dd, gsub("KEGG_", "", E1$kegg4$enrichRes@result$ID), ncol = 2,
	                    title = "Group 4", sub = "Negative control normalized",
	                    organism = organism, view_allpath = view_allpath, output=outputDir8)
	  if(!is.null(E1$kegg12$enrichRes) && nrow(E1$kegg12$enrichRes@result)>0)
	    arrangePathview(dd, gsub("KEGG_", "", E1$kegg12$enrichRes@result$ID), ncol = 2,
	                    title = "Group 1 & Group 2", sub = "Negative control normalized",
	                    organism=organism,view_allpath=view_allpath, output=outputDir8)
	  if(!is.null(E1$kegg13$enrichRes) && nrow(E1$kegg13$enrichRes@result)>0)
	    arrangePathview(dd, gsub("KEGG_", "", E1$kegg13$enrichRes@result$ID), ncol = 2,
	                    title = "Group 1 & Group 3", sub = "Negative control normalized",
	                    organism=organism, view_allpath=view_allpath, output=outputDir8)
	  if(!is.null(E1$kegg24$enrichRes) && nrow(E1$kegg24$enrichRes@result)>0)
	    arrangePathview(dd, gsub("KEGG_", "", E1$kegg24$enrichRes@result$ID), ncol = 2,
	                    title = "Group 2 & Group 4", sub = "Negative control normalized",
	                    organism=organism, view_allpath=view_allpath, output=outputDir8)
	  if(!is.null(E1$kegg34$enrichRes) && nrow(E1$kegg34$enrichRes@result)>0)
	    arrangePathview(dd, gsub("KEGG_", "", E1$kegg34$enrichRes@result$ID), ncol = 2,
	                    title = "Group 3 & Group 4", sub = "Negative control normalized",
	                    organism=organism, view_allpath=view_allpath, output=outputDir8)
  }
	dev.off()
}
