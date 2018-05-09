#' Do enrichment analysis using over-representation test
#'
#' Do enrichment analysis using over-representation test
#'
#' @docType methods
#' @name enrich.ORT
#' @rdname enrich.ORT
#' @aliases enrichORT
#'
#' @param gene A character vector, specifying the genelist to do enrichment analysis.
#' @param universe A character vector, specifying the backgound genelist, default is whole genome.
#' @param type Geneset category for testing, KEGG(default).
#' @param organism A character, specifying organism, such as "hsa" or "Human"(default), and "mmu" or "Mouse"
#' @param pvalueCutoff Pvalue cutoff.
#' @param qvalueCutoff Qvalue cutoff.
#' @param pAdjustMethod One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param minGSSize Minimal size of each geneSet for testing.
#' @param maxGSSize Maximal size of each geneSet for analyzing.
#'
#' @return A enrichResult instance.
#'
#' @author Wubing Zhang
#'
#' @note  See the vignette for an example of enrichment analysis using over-representation test
#' The source can be found by typing \code{MAGeCKFlute:::enrich.ORT}
#' or \code{getMethod("enrich.ORT")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/enrich.ORT.R}
#' Users should find it easy to customize this function.
#'
#' @seealso \code{\link{enrich.HGT}}
#' @seealso \code{\link{enrich.DAVID}}
#' @seealso \code{\link{enrich.GOstats}}
#' @seealso \code{\link{enrich.GSE}}
#' @seealso \code{\link{enrichment_analysis}}
#' @seealso \code{\link[DOSE]{enrichResult-class}}
#'
#' @examples
#' data(MLE_Data)
#' universe = TransGeneID(MLE_Data$Gene, "SYMBOL", "ENTREZID", organism = "hsa")
#' genes = universe[1:50]
#' enrichRes <- enrich.ORT(genes, universe)
#' head(enrichRes@result)
#'
#'
#' @import clusterProfiler
#' @export

enrich.ORT <- function(gene, universe=NULL, type="KEGG", organism = "hsa",pvalueCutoff = 1,
                       qvalueCutoff = 1, pAdjustMethod = "BH",minGSSize = 2, maxGSSize = 500){
  requireNamespace("clusterProfiler", quietly=TRUE) || stop("need clusterProfiler package")
  gene = unique(as.character(gene))
  if(!is.null(universe)) universe = unique(as.character(universe))
  organism = getOrg(organism)$org
  orgdb = getOrg(organism)$pkg
  #=======================
  if(type %in% c("BP", "CC", "MF")){
    enrichedRes = enrichGO(gene=gene,  universe=universe,  ont = type, OrgDb=orgdb,
                           pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff,
                           minGSSize=minGSSize, maxGSSize=maxGSSize)
    if(!is.null(enrichedRes))
        enrichedRes = simplify(enrichedRes, cutoff=0.7, by="p.adjust", select_fun=min)
  }
  if(type == "KEGG"){
    enrichedRes = enrichKEGG(gene=gene,  universe=universe, organism = organism,
                             pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff, qvalueCutoff=qvalueCutoff,
                             minGSSize = minGSSize, maxGSSize = maxGSSize)
  }
  if(type == "DO"){
    enrichedRes = enrichDO(gene=gene,  universe=universe, pAdjustMethod = pAdjustMethod,pvalueCutoff = pvalueCutoff,
                           qvalueCutoff=qvalueCutoff, minGSSize = minGSSize, maxGSSize = maxGSSize)
  }
  if(type == "MKEGG"){
    enrichedRes = enrichMKEGG(gene=gene,  universe=universe, organism = organism,pAdjustMethod=pAdjustMethod,
                              pvalueCutoff = pvalueCutoff,qvalueCutoff=qvalueCutoff, minGSSize = minGSSize, maxGSSize = maxGSSize)
  }
  if(type == "NCG"){
    enrichedRes = enrichNCG(gene=gene,  universe=universe, pAdjustMethod = pAdjustMethod,pvalueCutoff=pvalueCutoff,
                            qvalueCutoff=qvalueCutoff, minGSSize = minGSSize, maxGSSize = maxGSSize)
  }
  if(!is.null(enrichedRes)){
    geneID = strsplit(enrichedRes@result$geneID, "/")
    geneName = lapply(geneID, function(gid){
      SYMBOL = TransGeneID(gid, "ENTREZID", "SYMBOL", organism = organism)
      paste(SYMBOL, collapse = "/")
    })
    enrichedRes@result$geneName = unlist(geneName)
  }
  return(enrichedRes)
}

