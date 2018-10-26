#' Downstream analysis based on MAGeCK-MLE result
#'
#' Integrative analysis pipeline using the gene summary table in MAGeCK MLE results
#'
#' @docType methods
#' @name FluteMLE
#' @rdname FluteMLE
#' @aliases flutemle
#'
#' @param gene_summary Either a file path or a data frame, which contains columns of 'Gene',
#' \code{ctrlname}.beta and \code{treatname}.beta which corresponding to the parameter ctrlname and treatmname.
#' @param ctrlname A character vector, specifying the names of control samples.
#' @param treatname A character vector, specifying the names of treatment samples.
#' @param keytype Type of gene id in `gene_summary`, which should be one of "Entrez" or "Symbol".
#' @param organism A character, specifying organism, such as "hsa"("Human"), or "mmu"("Mouse").
#'
#' @param scale_cutoff Boolean or numeric, whether scale cutoff to whole genome level,
#' or how many standard deviation will be used as cutoff.
#' @param top An integer, specifying number of top selected genes to be labeled in rank figure.
#' @param bottom An integer, specifying number of bottom selected genes to be labeled in rank figure.
#' @param interestGenes A character vector, specifying interested genes to be labeled in rank figure.
#'
#' @param pathway_limit A two-length vector (default: c(3, 50)), specifying the min and
#' max size of pathways for enrichent analysis.
#' @param pvalueCutoff A numeric, specifying pvalue cutoff of enrichment analysis, default 1.
#' @param adjust One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param enrich_kegg One of "ORT"(Over-Representing Test), "GSEA"(Gene Set Enrichment Analysis), and "HGT"(HyperGemetric test).
#' @param gsea Boolean, indicating whether GSEA analysis is needed for positive and
#' negative selection genes.
#'
#' @param posControl A file path or a character vector, specifying a list of gene entrezid as
#' positive controls used for cell cycle normalization.
#' @param loess Boolean, whether include loess normalization in the pipeline.
#'
#' @param prefix A character, indicating the prefix of output file name, which can't contain special characters.
#' @param outdir Output directory on disk.
#' @param view_allpath Boolean, whether output all pathway view figures.
#'
#' @author Wubing Zhang
#'
#' @return All of the pipeline results is output into the \code{out.dir}/\code{prefix}_Results,
#' which includes a pdf file and many folders. The pdf file '{prefix}_Pipeline_results.pdf' is the
#' summary of pipeline results. For each section in this pipeline, figures and useful data are
#' outputed to corresponding subfolders.
#' Distribution_of_BetaScores: Density plot and violin plot of beta scores.
#' MAplot: Maplot for each normalized data.
#' Linear_Fitting_of_BetaScores: Linear fitting of beta scores indicates the difference of cell cycle
#' time between Control and Treatment samples.
#' Scatter_Treat_Ctrl: Positive selection and negative selection
#' Enrichment_Treat-Ctrl: Enrichment analysis for positive and negative selection genes
#' Pathview_Treat_Ctrl: Pathway view for top enriched pathways
#' Scatter_9Square: Using 9 Square to select drug related genes
#' Enrichment_9Square: Enrichment analysis for selected genes
#' Pathview_9Square: Pathway view for top enriched pathways
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
#' data(MLE_Data)
#' \dontrun{
#'   # functional analysis for MAGeCK MLE results
#'   FluteMLE(MLE_Data, ctrlname = c("D7_R1","D7_R2"), treatname = c("PLX7_R1","PLX7_R2"),
#'            prefix = "BRAF_D7", pvalueCutoff = 0.25, organism = "hsa")
#' }
#'
#' @import ggplot2 stats grDevices utils gridExtra grid
#' @export

FluteMLE <- function(gene_summary, ctrlname, treatname, keytype = "Symbol", organism = "hsa", # Input dataset
                     scale_cutoff = 1, top = 10, bottom = 10, interestGenes = NA, # Parameters for rank visualization
                     pathway_limit = c(3,50), pvalueCutoff=0.25, adjust="BH", enrich_kegg = "HGT", gsea = FALSE,
                     posControl = NULL, loess=FALSE,
                     prefix = "", outdir = ".", view_allpath = FALSE){

	#=========Prepare the running environment=========
	{
	  message(Sys.time(), " # Create output dir and pdf file...")
	  outdir = file.path(outdir, paste0(prefix, "Flute_Results"))
	  dir.create(file.path(outdir), showWarnings = FALSE)
    setwd(outdir)
	  output_pdf = paste0(prefix, "Flute.mle_summary.pdf")
	  if(loess){
	    pdf(output_pdf, width = 15, height = 7)
	  }else{
		  pdf(output_pdf, width=11, height = 6)
	  }
	  organism = getOrg(organism)$org
	}

	#=========Beta Score Preparation=========================
	{
	  beta = ReadBeta(gene_summary, keytype = "Symbol", organism = organism)
	  if(all(c(ctrlname, treatname) %in% colnames(beta))){
	    dd = beta[, c("Gene", "EntrezID", ctrlname, treatname)]
	  }else{
	    stop("No sample found!")
	  }
	  dd_essential = NormalizeBeta(dd, samples = c(ctrlname, treatname),
	                               method = "cell_cycle", posControl = posControl)
	  if(loess) dd_loess = NormalizeBeta(dd, samples = c(ctrlname, treatname), method = "loess")
	  rm(beta)
	}

	#========Distribution of all genes================================
	{
	  outputDir1 = "Distribution_of_BetaScores"
	  dir.create(outputDir1, showWarnings = FALSE)
	  idx_distr = c(ctrlname, treatname)
	  #Negative normalization
	  P1 = ViolinView(dd[, idx_distr], main = "Negative control normalized",
	                 filename = file.path(outputDir1, "violin_all_negative_normalized.png"))
	  P2 = DensityView(dd[, idx_distr], main = "Negative control normalized",
	                  filename = file.path(outputDir1, "density_all_negative_normalized.png"))
	  P3 = DensityDiffView(dd, ctrlname, treatname, main = "Negative control normalized",
	                  filename = file.path(outputDir1, "density_all_treat-ctrl_negative_normalized.png"))
	  #Essential gene normalization
	  P4 = ViolinView(dd_essential[,idx_distr], main = "Cell cycle normalized",
	                 filename = file.path(outputDir1, "violin_all_essential_normalized.png"))
	  P5 = DensityView(dd_essential[,idx_distr], main = "Cell cycle normalized",
	                  filename = file.path(outputDir1,"density_all_essential_normalized.png"))
	  P6 = DensityDiffView(dd_essential, ctrlname, treatname, main = "Cell cycle normalized",
	                  filename = file.path(outputDir1, "density_all_treat-ctrl_essential_normalized.png"))
	  #Loess normalization
	  if(loess){
  		P7 = ViolinView(dd_loess[, idx_distr], main = "Loess normalized",
  		               filename = file.path(outputDir1,"violin_all_loess_normalized.png"))
  		P8 = DensityView(dd_loess[, idx_distr],main = "Loess normalized",
  		                filename = file.path(outputDir1,"density_all_loess_normalized.png"))
  		P9 = DensityDiffView(dd_loess, ctrlname, treatname, main = "Loess normalized",
  		                filename = file.path(outputDir1, "density_all_treat-ctrl_loess_normalized.png"))
  		grid.arrange(P1, P4, P7, P2, P5, P8, ncol = 3)
	  }else grid.arrange(P1, P4, P2, P5, ncol = 2)
	  suppressWarnings(rm(P1, P2, P4, P5, P7, P8))
	}

	#========MAplot of treatment and control beta scores==================
	{
	  outputDir2 = "MAplot"
	  dir.create(outputDir2, showWarnings = FALSE)

	  P1 = MAView(dd, ctrlname, treatname, main = "Negative control normalized",
	          filename = file.path(outputDir2, "maplot_negative_normalized.png"))
	  P2 = MAView(dd_essential, ctrlname, treatname, main = "Cell cycle normalized",
	          filename = file.path(outputDir2, "maplot_essential_normalized.png"))
	  if(loess){
	    P4 = MAView(dd_loess, ctrlname, treatname, main = "Loess normalized",
	                    filename = file.path(outputDir2, "maplot_loess_normalized.png"))
	    grid.arrange(P3, P6, P9, P1, P2, P4, ncol = 3)
	  }else{
	    grid.arrange(P3, P6, P1, P2, ncol = 2)
	  }
	  suppressWarnings(rm(P1, P2, P3, P4, outputDir2))
	}

	#=============Distribution of essential genes====================
	{
	  outputDir3 = "Linear_Fitting_of_BetaScores"
	  dir.create(outputDir3, showWarnings = FALSE)
	  data(Zuber_Essential)
	  idx = which(rownames(dd) %in% Zuber_Essential$EntrezID)
	  #Negative control normalized
	  P1 = ViolinView(dd[idx, idx_distr], ylab = "Essential.B.S.", main = "Negative control normalized",
	                 filename = file.path(outputDir1, "violin_ess_negative_normalized.png"))
	  P2 = DensityView(dd[idx, idx_distr], xlab = "Essential.B.S.", main = "Negative control normalized",
	                  filename = file.path(outputDir1, "density_ess_negative_normalized.png"))
	  P3 = CellCycleView(dd[, idx_distr], ctrlname, treatname, main="Negative control normalized",
	                  filename = file.path(outputDir3, "Linear_all_negative_normalized.png"))
	  P4 = CellCycleView(dd[idx, idx_distr], ctrlname, treatname, main = "Negative control normalized",
	                  filename = file.path(outputDir3, "Linear_ess_negative_normalized.png"))
	  #Essential normalized
	  P5 = ViolinView(dd_essential[idx, idx_distr], ylab = "Essential.B.S.", main = "Cell cycle  normalized",
	                 filename = file.path(outputDir1, "violin_ess_essential_normalized.png"))
	  P6 = DensityView(dd_essential[idx, idx_distr], xlab = "Essential.B.S.", main = "Cell cycle  normalized",
	                  filename = file.path(outputDir1, "density_ess_essential_normalized.png"))
	  P7 = CellCycleView(dd_essential[, idx_distr], ctrlname, treatname, main = "Cell cycle  normalized",
	                  filename = file.path(outputDir3, "Linear_all_essential_normalized.png"))
	  P8 = CellCycleView(dd_essential[idx, idx_distr], ctrlname, treatname, main = "Cell cycle  normalized",
	                  filename = file.path(outputDir3, "Linear_ess_essential_normalized.png"))

	  #loess normalized
	  if(loess){
  		P9 = ViolinView(dd_loess[idx, idx_distr], ylab = "Essential.B.S.", main="Loess  normalized",
  		               filename = file.path(outputDir1, "violin_ess_loess_normalized.png"))
  		P10 = DensityView(dd_loess[idx,idx_distr], xlab = "Essential.B.S.", main = "Loess  normalized",
  		                 filename = file.path(outputDir1, "density_ess_loess_normalized.png"))
  		P11 = CellCycleView(dd_loess[,idx_distr], ctrlname, treatname, main = "Loess  normalized",
  		                 filename = file.path(outputDir3, "Linear_all_loess_normalized.png"))
  		P12 = CellCycleView(dd_loess[idx,idx_distr], ctrlname, treatname, main = "Loess  normalized",
  		                 filename = file.path(outputDir3, "Linear_ess_loess_normalized.png"))

  		grid.arrange(P1, P5, P9, P2, P6, P10, ncol = 3)
  		grid.arrange(P3, P7, P11, P4, P8, P12, ncol = 3)
	  }else{
  		grid.arrange(P1, P5, P2, P6, ncol = 2)
  		grid.arrange(P3, P7, P4, P8, ncol = 2)
	  }
	  suppressWarnings(rm(P1, P2, P3, P4, P5, P6, P7, P8, P9, P10, P11, P12, outputDir3, idx))
	}

  dd$Control = rowMeans(dd[, ctrlname, drop = FALSE])
  dd$Treatment = rowMeans(dd[, treatname, drop = FALSE])
  dd.diff = dd$Treatment - dd$Control
  names(dd.diff) = dd$Gene
  dd_essential$Control = rowMeans(dd_essential[, ctrlname, drop = FALSE])
  dd_essential$Treatment = rowMeans(dd_essential[, treatname, drop = FALSE])
  dd_essential.diff = dd_essential$Treatment - dd_essential$Control
  names(dd_essential.diff) = dd_essential$Gene
  if(loess){
    dd_loess$Control = rowMeans(dd_loess[, ctrlname, drop = FALSE])
    dd_loess$Treatment = rowMeans(dd_loess[, treatname, drop = FALSE])
    dd_loess.diff = dd_loess$Treatment - dd_loess$Control
    names(dd_loess.diff) = dd_loess$Gene
  }

	# =========Drug-targeted genes=================================
	{
	  outputDir4 = "Scatter_Treat_Ctrl"
	  dir.create(outputDir4, showWarnings = FALSE)

	  # Negative normalized
	  P1 = ScatterView(dd, main = "Negative control normalized", scale_cutoff = scale_cutoff,
	                filename = file.path(outputDir4,"Scatter_Treat-Ctrl_negative_normalized.png"))
	  P2 = RankView(dd.diff, genelist = interestGenes, top = top, bottom = bottom, main = "Negative control normalized",
	              cutoff = c(-CutoffCalling(dd.diff, scale = scale_cutoff), CutoffCalling(dd.diff, scale = scale_cutoff)),
	              filename = file.path(outputDir4, "Rank_Treat-Ctrl_negative_normalized.png"))
	  #Essential normalized
	  P3=ScatterView(dd_essential, main="Cell cycle normalized", scale_cutoff = scale_cutoff,
	               filename = file.path(outputDir4, "Scatter_Treat-Ctrl_essential_normalized.png"))
	  P4=RankView(dd_essential.diff, genelist = interestGenes, top = top, bottom = bottom, main = "Cell cycle  normalized",
	              cutoff = c(-CutoffCalling(dd_essential.diff, scale = scale_cutoff),
	                     CutoffCalling(dd_essential.diff, scale = scale_cutoff)),
	              filename = file.path(outputDir4, "Rank_Treat-Ctrl_essential_normalized.png"))

	  # Loess normalized
	  if(loess){
  		P5 = ScatterView(dd_loess, main = "Loess normalized",
  		             scale_cutoff = scale_cutoff,
  		             filename = file.path(outputDir4, "Scatter_Treat-Ctrl_loess_normalized.png"))
  		P6 = RankView(dd_loess.diff, genelist = interestGenes, top = top, bottom = bottom, main = "Loess  normalized",
  		            cutoff = c(-CutoffCalling(dd_loess.diff, scale = scale_cutoff),
  		                     CutoffCalling(dd_loess.diff, scale = scale_cutoff)),
  		          filename = file.path(outputDir4, "Rank_Treat-Ctrl_loess_normalized.png"))

  		grid.arrange(P1, P3, P5, P2, P4, P6, ncol = 3)
	  }else{
	    grid.arrange(P1, P3, P2, P4, ncol = 2)
	  }
	}

	#==========Enrichment AB group genes=========================
	{
	  outputDir5 = "Enrichment_Treat-Ctrl/"
	  dir.create(outputDir5, showWarnings=FALSE)

	  E1 = EnrichAB(P1$data, pvalue = pvalueCutoff, enrich_method = enrich_kegg,
	                organism = organism, adjust = adjust, gsea = gsea,
	                filename = "Negative_ctrl_normalized", out.dir = outputDir5)
	  E2 = EnrichAB(P3$data, pvalue = pvalueCutoff,
	                enrich_method = enrich_kegg, organism=organism,
	                adjust=adjust, gsea=gsea,
	                filename="Essential_normalized", out.dir=outputDir5)

	  if(loess){
  		E3 = EnrichAB(P5$data, pvalue=pvalueCutoff,
  		              enrich_method = enrich_kegg, organism=organism,
  		              adjust=adjust, gsea=gsea,filename="Loess_normalized",
  		              out.dir=outputDir5)

  		# grid.arrange(E1$keggA$gridPlot, E2$keggA$gridPlot, E3$keggA$gridPlot,
  		#              E1$keggB$gridPlot, E2$keggB$gridPlot, E3$keggB$gridPlot, ncol = 3)
  		grid.arrange(E1$keggA$gridPlot, E2$keggA$gridPlot, E3$keggA$gridPlot,
  		             E1$bpA$gridPlot, E2$bpA$gridPlot, E3$bpA$gridPlot,
  		             ncol = 3)
  		if(gsea) grid.arrange(E1$gseA$gridPlot, E2$gseA$gridPlot,
  		                      E3$gseA$gridPlot, E1$gseA$gseaplot,
  		                      E2$gseA$gseaplot, E3$gseA$gseaplot, ncol = 3)
  		grid.arrange(E1$keggB$gridPlot, E2$keggB$gridPlot, E3$keggB$gridPlot,
  		             E1$bpB$gridPlot, E2$bpB$gridPlot,
  		             E3$bpB$gridPlot, ncol = 3)
  		if(gsea)grid.arrange(E1$gseB$gridPlot, E2$gseB$gridPlot,
  		                     E3$gseB$gridPlot, E1$gseB$gseaplot,
  		                     E2$gseB$gseaplot, E3$gseB$gseaplot, ncol = 3)
	  }else{
      # grid.arrange(E1$keggA$gridPlot, E2$keggA$gridPlot, E1$keggB$gridPlot, E2$keggB$gridPlot, ncol = 2)
	    grid.arrange(E1$keggA$gridPlot, E2$keggA$gridPlot, E1$bpA$gridPlot,
	                 E2$bpA$gridPlot, ncol = 2)
  		if(gsea) grid.arrange(E1$gseA$gridPlot, E2$gseA$gridPlot,
  		                      E1$gseA$gseaplot, E2$gseA$gseaplot, ncol = 2)
  		grid.arrange(E1$keggB$gridPlot, E2$keggB$gridPlot, E1$bpB$gridPlot,
  		             E2$bpB$gridPlot, ncol = 2)
  		if(gsea) grid.arrange(E1$gseB$gridPlot, E2$gseB$gridPlot,
  		                      E1$gseB$gseaplot, E2$gseB$gseaplot, ncol = 2)
	  }
	}

	#======Pathview for GroupA and GroupB significant enriched pathways
	{
	  dir.create("Pathview_Treat_Ctrl", showWarnings=FALSE)
	  #Pathway view for top 4 pathway
	  if(!is.null(E1$keggA$enrichRes) && nrow(E1$keggA$enrichRes@result)>0)
  	  arrangePathview(dd, E1$keggA$enrichRes@result$ID, "Group A",
  	              "Negative control normalized",organism=organism,
  	              view_allpath=view_allpath, output="Pathview_Treat_Ctrl")
	  if(!is.null(E1$keggB$enrichRes) && nrow(E1$keggB$enrichRes@result)>0)
  	  arrangePathview(dd, E1$keggB$enrichRes@result$ID, "Group B",
  	              "Negative control normalized",organism=organism,
  	              view_allpath=view_allpath, output="Pathview_Treat_Ctrl")
	  if(!is.null(E2$keggA$enrichRes) && nrow(E2$keggA$enrichRes@result)>0)
  	  arrangePathview(dd_essential, E2$keggA$enrichRes@result$ID, "Group A",
  	              "Cell cycle normalized",organism=organism,
  	              view_allpath=view_allpath, output="Pathview_Treat_Ctrl")
	  if(!is.null(E2$keggB$enrichRes) && nrow(E2$keggB$enrichRes@result)>0)
  	  arrangePathview(dd_essential, E2$keggB$enrichRes@result$ID, "Group B",
  	              "Cell cycle normalized",organism=organism,
  	              view_allpath=view_allpath, output="Pathview_Treat_Ctrl")

	  if(loess){
	    if(!is.null(E3$keggA$enrichRes) && nrow(E3$keggA$enrichRes@result)>0)
    		arrangePathview(dd_loess, E3$keggA$enrichRes@result$ID, "Group A",
    		            "Loess normalized",organism=organism,
    		            view_allpath=view_allpath, output="Pathview_Treat_Ctrl")
	    if(!is.null(E3$keggB$enrichRes) && nrow(E3$keggB$enrichRes@result)>0)
    		arrangePathview(dd_loess, E3$keggB$enrichRes@result$ID, "Group B",
    		            "Loess normalized",organism=organism,
    		            view_allpath=view_allpath, output="Pathview_Treat_Ctrl")
	  }
	  suppressWarnings(rm(E1, E2, E3, outputDir5))
	}

	# ===============9 squares=====================================
	{
	  dir.create("Scatter_9Square", showWarnings=FALSE)
	  outputDir6 <- "Scatter_9Square/Square9_"

	  P1 = SquareView(dd, label="Gene", main="Negative control normalized", scale_cutoff=scale_cutoff,
	                 filename=paste0(outputDir6,"scatter_negative_normalized.png"))

	  P2 = SquareView(dd_essential, label="Gene", main="Cell cycle normalized", scale_cutoff=scale_cutoff,
	                filename=paste0(outputDir6,"scatter_cellcycle_normalized.png"))

	  if(loess){
  		P3 = SquareView(dd_loess, label="Gene", main="Loess normalized", scale_cutoff=scale_cutoff,
  		              filename=paste0(outputDir6,"scatter_loess_normalized.png"))
  		grid.arrange(P1,P2,P3, ncol = 3)
	  }else{grid.arrange(P1,P2, ncol = 2, bottom="", left="", right="", top="")}
	}

	#==============9 Square grouped gene enrichment======================
	{
	  dir.create("Enrichment_9Square", showWarnings=FALSE)
	  E1 = EnrichSquare(P1$data, pvalue = pvalueCutoff,adjust = adjust,
	                    enrich_method = enrich_kegg, organism=organism,
        						 filename="negative_normalized",
        						 out.dir="Enrichment_9Square")
	  E2 = EnrichSquare(P2$data, pvalue=pvalueCutoff,adjust=adjust, enrich_method=enrich_kegg,
	                    organism=organism, filename="essential_normalized",
        						 out.dir = "Enrichment_9Square")

	  if(loess){
  		E3 = EnrichSquare(P3$data, pvalue=pvalueCutoff,adjust=adjust, organism=organism,
  		                  enrich_method=enrich_kegg, filename="loess_normalized",
        						   out.dir="Enrichment_9Square")

  		grid.arrange(E1$kegg1$gridPlot, E2$kegg1$gridPlot, E3$kegg1$gridPlot,
  		             E1$bp1$gridPlot, E2$bp1$gridPlot, E3$bp1$gridPlot, ncol = 3)
  		grid.arrange(E1$kegg2$gridPlot, E2$kegg2$gridPlot, E3$kegg2$gridPlot,
  		             E1$bp2$gridPlot, E2$bp2$gridPlot, E3$bp2$gridPlot, ncol = 3)
  		grid.arrange(E1$kegg3$gridPlot, E2$kegg3$gridPlot, E3$kegg3$gridPlot,
  		             E1$bp3$gridPlot, E2$bp3$gridPlot, E3$bp3$gridPlot, ncol = 3)
  		grid.arrange(E1$kegg4$gridPlot, E2$kegg4$gridPlot, E3$kegg4$gridPlot,
  		             E1$bp4$gridPlot, E2$bp4$gridPlot, E3$bp4$gridPlot, ncol = 3)
  		grid.arrange(E1$kegg13$gridPlot, E2$kegg13$gridPlot, E3$kegg13$gridPlot,
  		             E1$bp13$gridPlot, E2$bp13$gridPlot, E3$bp13$gridPlot, ncol = 3)
  		grid.arrange(E1$kegg14$gridPlot, E2$kegg14$gridPlot, E3$kegg14$gridPlot,
  		             E1$bp14$gridPlot, E2$bp14$gridPlot, E3$bp14$gridPlot, ncol = 3)
  		grid.arrange(E1$kegg23$gridPlot, E2$kegg23$gridPlot, E3$kegg23$gridPlot,
  		             E1$bp23$gridPlot, E2$bp23$gridPlot, E3$bp23$gridPlot, ncol = 3)
  		grid.arrange(E1$kegg24$gridPlot, E2$kegg24$gridPlot, E3$kegg24$gridPlot,
  		             E1$bp24$gridPlot, E2$bp24$gridPlot, E3$bp24$gridPlot, ncol = 3)
  		}else{
  		grid.arrange(E1$kegg1$gridPlot, E2$kegg1$gridPlot, E1$bp1$gridPlot,
  		             E2$bp1$gridPlot, ncol = 2)
  		grid.arrange(E1$kegg2$gridPlot, E2$kegg2$gridPlot, E1$bp2$gridPlot,
  		             E2$bp2$gridPlot, ncol = 2)
  		grid.arrange(E1$kegg3$gridPlot, E2$kegg3$gridPlot, E1$bp3$gridPlot,
  		             E2$bp3$gridPlot, ncol = 2)
  		grid.arrange(E1$kegg4$gridPlot, E2$kegg4$gridPlot, E1$bp4$gridPlot,
  		             E2$bp4$gridPlot, ncol = 2)
  		grid.arrange(E1$kegg13$gridPlot, E2$kegg13$gridPlot, E1$bp13$gridPlot,
  		             E2$bp13$gridPlot, ncol = 2)
  		grid.arrange(E1$kegg14$gridPlot, E2$kegg14$gridPlot, E1$bp14$gridPlot,
  		             E2$bp14$gridPlot, ncol = 2)
  		grid.arrange(E1$kegg23$gridPlot, E2$kegg23$gridPlot, E1$bp23$gridPlot,
  		             E2$bp23$gridPlot, ncol = 2)
  		grid.arrange(E1$kegg24$gridPlot, E2$kegg24$gridPlot, E1$bp24$gridPlot,
  		             E2$bp24$gridPlot, ncol = 2)
		}
	}

	#========Pathway view for 9 square genesets================
	{
	  dir.create("Pathview_9Square", showWarnings=FALSE)
	  if(!is.null(E1$kegg1$enrichRes) && nrow(E1$kegg1$enrichRes@result)>0)
  	  arrangePathview(dd, E1$kegg1$enrichRes@result$ID,"Group 1",
  	              "Negative control normalized",organism=organism, view_allpath=view_allpath,
  	              output="Pathview_9Square")
	  if(!is.null(E1$kegg2$enrichRes) && nrow(E1$kegg2$enrichRes@result)>0)
  	  arrangePathview(dd, E1$kegg2$enrichRes@result$ID,"Group 2",
  	              "Negative control normalized",organism=organism,view_allpath=view_allpath,
  	              output="Pathview_9Square")
	  if(!is.null(E1$kegg3$enrichRes) && nrow(E1$kegg3$enrichRes@result)>0)
  	  arrangePathview(dd, E1$kegg3$enrichRes@result$ID,"Group 3",
  	              "Negative control normalized",organism=organism,view_allpath=view_allpath,
  	              output="Pathview_9Square")
	  if(!is.null(E1$kegg4$enrichRes) && nrow(E1$kegg4$enrichRes@result)>0)
  	  arrangePathview(dd, E1$kegg4$enrichRes@result$ID,"Group 4",
  	              "Negative control normalized",organism=organism,view_allpath=view_allpath,
  	              output="Pathview_9Square")
	  if(!is.null(E1$kegg13$enrichRes) && nrow(E1$kegg13$enrichRes@result)>0)
  	  arrangePathview(dd, E1$kegg13$enrichRes@result$ID,"Group 1 & Group 3",
  	              "Negative control normalized",organism=organism,view_allpath=view_allpath,
  	              output="Pathview_9Square")
	  if(!is.null(E1$kegg14$enrichRes) && nrow(E1$kegg14$enrichRes@result)>0)
  	  arrangePathview(dd, E1$kegg14$enrichRes@result$ID,"Group 1 & Group 4",
  	              organism=organism, view_allpath=view_allpath,
  	              output="Pathview_9Square")
	  if(!is.null(E1$kegg23$enrichRes) && nrow(E1$kegg23$enrichRes@result)>0)
  	  arrangePathview(dd, E1$kegg23$enrichRes@result$ID,"Group 2 & Group 3",
  	              "Negative control normalized",organism=organism,
  	              view_allpath=view_allpath, output="Pathview_9Square")
	  if(!is.null(E1$kegg24$enrichRes) && nrow(E1$kegg24$enrichRes@result)>0)
  	  arrangePathview(dd, E1$kegg24$enrichRes@result$ID,
  	              "Group 2 & Group 4","Negative control normalized",
  	              organism=organism, view_allpath=view_allpath,
  	              output="Pathview_9Square")
	  if(!is.null(E2$kegg1$enrichRes) && nrow(E2$kegg1$enrichRes@result)>0)
	    arrangePathview(dd_essential, E2$kegg1$enrichRes@result$ID,"Group 1",
	                    "Cell cycle normalized",organism=organism,view_allpath=view_allpath,
	                    output="Pathview_9Square")
	  if(!is.null(E2$kegg2$enrichRes) && nrow(E2$kegg2$enrichRes@result)>0)
	    arrangePathview(dd_essential, E2$kegg2$enrichRes@result$ID,"Group 2",
	                    "Cell cycle normalized",organism=organism,view_allpath=view_allpath,
	                    output="Pathview_9Square")
	  if(!is.null(E2$kegg3$enrichRes) && nrow(E2$kegg3$enrichRes@result)>0)
	    arrangePathview(dd_essential, E2$kegg3$enrichRes@result$ID,"Group 3",
	                    "Cell cycle normalized",organism=organism,view_allpath=view_allpath,
	                    output="Pathview_9Square")
	  if(!is.null(E2$kegg4$enrichRes) && nrow(E2$kegg4$enrichRes@result)>0)
	    arrangePathview(dd_essential, E2$kegg4$enrichRes@result$ID,"Group 4",
	                    "Cell cycle normalized",organism=organism,view_allpath=view_allpath,
	                    output="Pathview_9Square")
	  if(!is.null(E2$kegg13$enrichRes) && nrow(E2$kegg13$enrichRes@result)>0)
	    arrangePathview(dd_essential, E2$kegg13$enrichRes@result$ID,
	                    "Group 1 & Group 3","Cell cycle normalized",
	                    organism=organism, view_allpath=view_allpath,
	                    output="Pathview_9Square")
	  if(!is.null(E2$kegg14$enrichRes) && nrow(E2$kegg14$enrichRes@result)>0)
	    arrangePathview(dd_essential, E2$kegg14$enrichRes@result$ID,
	                    "Group 1 & Group 4","Cell cycle normalized",
	                    organism=organism, view_allpath=view_allpath,
	                    output="Pathview_9Square")
	  if(!is.null(E2$kegg23$enrichRes) && nrow(E2$kegg23$enrichRes@result)>0)
	    arrangePathview(dd_essential, E2$kegg23$enrichRes@result$ID,
	                    "Group 2 & Group 3","Cell cycle normalized",
	                    organism=organism, view_allpath=view_allpath,
	                    output="Pathview_9Square")
	  if(!is.null(E2$kegg24$enrichRes) && nrow(E2$kegg24$enrichRes@result)>0)
	    arrangePathview(dd_essential, E2$kegg24$enrichRes@result$ID,
	                    "Group 2 & Group 4","Cell cycle normalized",
	                    organism=organism, view_allpath=view_allpath,
	                    output="Pathview_9Square")

	  if(loess){
	    if(!is.null(E3$kegg1$enrichRes) && nrow(E3$kegg1$enrichRes@result)>0)
    		arrangePathview(dd_loess, E3$kegg1$enrichRes@result$ID,"Group 1",
    		            "Loess normalized",organism=organism,
    		            view_allpath=view_allpath, output="Pathview_9Square")
	    if(!is.null(E3$kegg2$enrichRes) && nrow(E3$kegg2$enrichRes@result)>0)
    		arrangePathview(dd_loess, E3$kegg2$enrichRes@result$ID,"Group 2",
    		            "Loess normalized",organism=organism,
    		            view_allpath=view_allpath, output="Pathview_9Square")
	    if(!is.null(E3$kegg3$enrichRes) && nrow(E3$kegg3$enrichRes@result)>0)
    		arrangePathview(dd_loess, E3$kegg3$enrichRes@result$ID,"Group 3",
    		            "Loess normalized",organism=organism,
    		            view_allpath=view_allpath, output="Pathview_9Square")
	    if(!is.null(E3$kegg4$enrichRes) && nrow(E3$kegg4$enrichRes@result)>0)
    		arrangePathview(dd_loess, E3$kegg4$enrichRes@result$ID,"Group 4",
    		            "Loess normalized",organism=organism,
    		            view_allpath=view_allpath, output="Pathview_9Square")
	    if(!is.null(E3$kegg13$enrichRes) && nrow(E3$kegg13$enrichRes@result)>0)
	      arrangePathview(dd_loess, E3$kegg13$enrichRes@result$ID,
    		            "Group 1 & Group 3","Loess normalized",organism=organism,
    		            view_allpath=view_allpath, output="Pathview_9Square")
	    if(!is.null(E3$kegg14$enrichRes) && nrow(E3$kegg14$enrichRes@result)>0)
    		arrangePathview(dd_loess, E3$kegg14$enrichRes@result$ID,
    		            "Group 1 & Group 4","Loess normalized",organism=organism,
    		            view_allpath=view_allpath, output="Pathview_9Square")
	    if(!is.null(E3$kegg23$enrichRes) && nrow(E3$kegg23$enrichRes@result)>0)
  	    arrangePathview(dd_loess, E3$kegg23$enrichRes@result$ID,
      		            "Group 2 & Group 3","Loess normalized",organism=organism,
      		            view_allpath=view_allpath, output="Pathview_9Square")
	    if(!is.null(E3$kegg24$enrichRes) && nrow(E3$kegg24$enrichRes@result)>0)
  	    arrangePathview(dd_loess, E3$kegg24$enrichRes@result$ID,
      		            "Group 2 & Group 4","Loess normalized",organism=organism,
      		            view_allpath=view_allpath, output="Pathview_9Square")
	  }
	}
	dev.off()
}
