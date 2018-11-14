#' Enrichment analysis
#'
#' Enrichment analysis
#'
#' @docType methods
#' @name enrichment_analysis
#' @rdname enrichment_analysis
#' @aliases enrichment
#'
#' @param geneList A numeric vector with gene as names.
#' @param universe A character vector, specifying the backgound genelist, default is whole genome.
#' @param method One of "ORT"(Over-Representing Test), "GSEA"(Gene Set Enrichment Analysis), and "HGT"(HyperGemetric test).
#' @param keytype "Entrez" or "Symbol".
#' @param type Geneset category for testing, one of 'GOBP+GOMF' (default), 'GOBP', 'GOMF', 'GOCC',
#' 'KEGG', 'BIOCARTA', 'REACTOME', 'TFT', 'IMMUNOLOGIC', 'ONCOGENIC',
#' or 'All' and any combination of them, such as 'IMMUNOLOGIC+GOBP+KEGG'.
#' @param organism 'hsa' or 'mmu'.
#' @param pvalueCutoff Pvalue cutoff.
#' @param pAdjustMethod One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param limit A two-length vector (default: c(3, 50)), specifying the minimal and
#' maximal size of gene sets for enrichent analysis.
#' @param plotTitle Same as 'title' in 'plot'.
#' @param color Color of points.
#'
#' @return A list, including two items, \code{gridPlot} and \code{enrichRes}. \code{gridPlot} is
#' a ggplot object, and \code{enrichRes} is a enrichResult instance.
#'
#' @author Feizhen Wu
#'
#' @seealso \code{\link{enrich.GSE}}
#' @seealso \code{\link{enrich.ORT}}
#' @seealso \code{\link{enrich.HGT}}
#' @seealso \code{\link[DOSE]{enrichResult-class}}
#'
#' @examples
#' data(geneList, package = "DOSE")
#' genes <- geneList[1:100]
#' keggA = enrichment_analysis(genes, method = "HGT", type = "KEGG")
#' print(keggA$gridPlot)
#'
#' @import DOSE
#' @export

enrichment_analysis = function(geneList, universe = NULL, method = "ORT", keytype = "Entrez",
                               type = "KEGG", organism = "hsa", pvalueCutoff = 0.25, pAdjustMethod = "BH",
                               limit = c(3, 50), plotTitle = NULL, color = "#3f90f7"){

  requireNamespace("stats", quietly=TRUE) || stop("need stats package")
  result = list()
  type = toupper(type[1])
  methods = c("ORT", "GSEA", "HGT")
  names(methods) = toupper(methods)
  if(class(method)=="character") method = toupper(method)
  method = methods[method]

  # The number of input genes < limit
  if(length(geneList) < limit[1] + 1){
    p1 = ggplot()
    p1 = p1 + geom_text(aes(x=0, y=0, label = "Too less input genes"), size = 6)
    p1 = p1 + labs(title = plotTitle)
    p1 = p1 + theme(plot.title = element_text(size=10))
    p1 = p1 + theme_void()
    p1 = p1 + theme(plot.title = element_text(hjust = 0.5))
    result$enrichRes = NULL
    result$gridPlot = p1
    return(result)
  }

  # Gene Set Enrichment Analysis
  if(method == "GSEA"){
    message(Sys.time(), " # Running GSEA for ", type)
    geneList = geneList[order(geneList,decreasing = TRUE)]
    enrichRes <- enrich.GSE(geneList, type = type, pvalueCutoff = pvalueCutoff,
                            pAdjustMethod = pAdjustMethod, organism = organism, limit = limit)
  }

  ## Over-Representation Analysis
  if(method == "ORT"){
    message(Sys.time(), " # Running Over-Representation Test for ", type)
    enrichRes <- enrich.ORT(geneList, universe = universe, keytype = keytype, type = type,
                            organism = organism, pvalueCutoff = pvalueCutoff,
                            pAdjustMethod = pAdjustMethod, limit = limit)
  }

  ## HyperGeometric test
  if(method == "HGT"){
    message(Sys.time(), " # Running Hypergeometric test for ", type)
    enrichRes = enrich.HGT(geneList, universe, keytype = keytype, type = type, organism = organism,
                           pvalueCutoff  = pvalueCutoff, pAdjustMethod = pAdjustMethod, limit = limit)
  }

  ## Visualization of enriched pathways
  result$enrichRes = enrichRes
  if(!is.null(enrichRes)){
    gridPlot <- EnrichedGSEView(enrichRes@result, plotTitle = plotTitle, color = color)
  }else{ gridPlot = noEnrichPlot(plotTitle) }
  result$gridPlot = gridPlot

  return(result)
}

