#' Enrichment analysis
#'
#' Enrichment analysis
#'
#' @docType methods
#' @name enrichment_analysis
#' @rdname enrichment_analysis
#' @aliases enrichment
#'
#' @param geneList A character vector or a ranked numeric vector(for GSEA) with names of geneid,
#' specifying the genelist to do enrichment analysis.
#' @param universe A character vector, specifying the backgound genelist, default is whole genome.
#' @param method One of "ORT"(Over-Representing Test), "GSEA"(Gene Set Enrichment Analysis), "DAVID",
#' "GOstats", and "HGT"(HyperGemetric test), or index from 1 to 5
#' @param type Geneset category for testing, KEGG(default).
#' @param organism A character, specifying organism, such as "hsa" or "Human"(default), and "mmu" or "Mouse"
#' @param pvalueCutoff Pvalue cutoff.
#' @param qvalueCutoff Qvalue cutoff.
#' @param pAdjustMethod One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param minGSSize Minimal size of each geneSet for testing.
#' @param maxGSSize Maximal size of each geneSet for analyzing.
#' @param plotTitle Same as 'title' in 'plot'.
#' @param color Color of points.
#'
#' @return A list, including two items, \code{gridPlot} and \code{enrichRes}. \code{gridPlot} is
#' a ggplot object, and \code{enrichRes} is a enrichResult instance.
#'
#' @author Feizhen Wu
#'
#' @note  See the vignette for an example of enrichment analysis
#' The source can be found by typing \code{MAGeCKFlute:::enrichment_analysis}
#' or \code{getMethod("enrichment_analysis")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/enrichment_analysis.R}
#' Users should find it easy to customize this function.
#'
#' @seealso \code{\link{enrich.GOstats}}
#' @seealso \code{\link{enrich.DAVID}}
#' @seealso \code{\link{enrich.GSE}}
#' @seealso \code{\link{enrich.ORT}}
#' @seealso \code{\link{enrich.HGT}}
#' @seealso \code{\link[DOSE]{enrichResult-class}}
#'
#' @examples
#' data(MLE_Data)
#' universe = TransGeneID(MLE_Data$Gene, "SYMBOL", "ENTREZID", organism = "hsa")
#' genes = universe[1:50]
#' keggA = enrichment_analysis(geneList = genes, universe=universe, method = "HGT",
#'                           type = "KEGG", organism = "hsa", color="#6daf61")
#' print(keggA$gridPlot)
#'
#'
#' @export

#====enrichment analysis===================================
enrichment_analysis = function(geneList, universe=NULL, method=1, type="KEGG", organism="hsa",
                               pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
                               minGSSize = 2, maxGSSize = 500, plotTitle=NULL, color="#3f90f7"){

  requireNamespace("stats", quietly=TRUE) || stop("need stats package")
  result=list()
  type=toupper(type[1])
  methods = c("ORT", "GSEA", "DAVID", "GOstats", "HGT")
  names(methods) = toupper(methods)
  if(class(method)=="character"){method = toupper(method)}
  method = methods[method]
  #======================================================================================
  if(length(geneList)<10){
    p1=ggplot()
    p1=p1+geom_text(aes(x=0,y=0,label="Less than 10 genes"),size=6)
    p1=p1+labs(title=plotTitle)
    p1=p1+theme(plot.title = element_text(size=10))
    p1=p1+theme_void()
    p1=p1+theme(plot.title = element_text(hjust = 0.5))

    result$enrichRes = NULL
    result$gridPlot = p1
    return(result)
  }
  #====Gene Set Enrichment Analysis=======================================================
  if(method == "GSEA"){
    loginfo(paste0('Running GSEA for ', type))
    geneList=geneList[order(geneList,decreasing = TRUE)]
    enrichRes <- enrich.GSE(geneList, type = type, pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod,
                            organism=organism, minGSSize = minGSSize, maxGSSize = maxGSSize)
    result$enrichRes = enrichRes
    if(!is.null(enrichRes)){
      gridPlot <- EnrichedGSEView(enrichRes@result, plotTitle, color=color)
    }else{
      p1=ggplot()
      p1=p1+geom_text(aes(x=0,y=0,label="No enriched terms"),size=6)
      p1=p1+labs(title=plotTitle)
      p1=p1+theme(plot.title = element_text(size=10))
      p1=p1+theme_void()
      p1=p1+theme(plot.title = element_text(hjust = 0.5))
      gridPlot = p1
    }
    result$gridPlot = gridPlot
    return(result)
  }
  #====Over-Representation Analysis======================================================
  if(method == "ORT"){
    loginfo(paste0('Running Over-Representation Test for ', type))
    enrichRes <- enrich.ORT(gene = geneList, universe = universe, type = type, organism=organism,
                            pvalueCutoff=pvalueCutoff, qvalueCutoff = qvalueCutoff, pAdjustMethod = pAdjustMethod,
                            minGSSize = minGSSize, maxGSSize = maxGSSize)
  }
  #=============DAVID=======================================================================
  if(method == "DAVID"){
    loginfo(paste0('Running DAVID for ', type))
    if(type == "KEGG"){
      annotation = "KEGG_PATHWAY"
    }else if(type %in% c("BP", "CC", "MF")){
      annotation = paste("GOTERM", type, "FAT", sep="_")
    }else if(type == "DO"){
      annotation = "OMIM_DISEASE"
    }else{annotation = type}

    enrichRes = enrich.DAVID(gene = geneList, universe = universe,
                             minGSSize = minGSSize, maxGSSize = maxGSSize, annotation  = annotation,
                             pvalueCutoff  = pvalueCutoff, pAdjustMethod = pAdjustMethod, qvalueCutoff= qvalueCutoff)
  }
  #===============================GOstats enrichment analysis============================
  if(method == "GOstats"){
    loginfo(paste0('Running GOstats test for ', type))
    enrichRes = enrich.GOstats(gene = geneList, universe = universe, type  = type, organism=organism,
                             pvalueCutoff  = pvalueCutoff, pAdjustMethod = pAdjustMethod)
  }
  #==================================HyperGeometric test=================================
  if(method == "HGT"){
    loginfo(paste0('Running Hypergeometric test for ', type))
    enrichRes = enrich.HGT(gene = geneList, universe = universe, type  = type, organism=organism,
                             pvalueCutoff  = pvalueCutoff, pAdjustMethod = pAdjustMethod,
                           minGSSize = minGSSize, maxGSSize = maxGSSize)
  }

  result$enrichRes = enrichRes
  if(!is.null(enrichRes)){
    gridPlot <- EnrichedView(enrichRes@result, plotTitle, color=color)
  }else{
    p1=ggplot()
    p1=p1+geom_text(aes(x=0,y=0,label="No enriched terms"),size=6)
    p1=p1+labs(title=plotTitle)
    p1=p1+theme(plot.title = element_text(size=10))
    p1=p1+theme_void()
    p1=p1+theme(plot.title = element_text(hjust = 0.5))
    gridPlot = p1
  }

  result$gridPlot = gridPlot

  return(result)
}

