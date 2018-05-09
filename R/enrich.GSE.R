#' GSEA
#'
#' A universal gene set enrichment analysis tools
#'
#' @docType methods
#' @name enrich.GSE
#' @rdname enrich.GSE
#' @aliases enrichGSE
#'
#' @param geneList A order ranked numeric vector with geneid as names.
#' @param type A character, indicating geneset category for testing, "KEGG"(default).
#' @param organism A character, specifying organism, such as "hsa" or "Human"(default), and "mmu" or "Mouse"
#' @param minGSSize Minimal size of each geneSet for testing.
#' @param maxGSSize Maximal size of each geneSet for analyzing.
#' @param pvalueCutoff Pvalue cutoff.
#' @param pAdjustMethod One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#'
#' @return A enrichResult instance.
#'
#' @author Wubing Zhang
#'
#' @note  See the vignette for an example of GSEA
#' The source can be found by typing \code{MAGeCKFlute:::enrich.GSE}
#' or \code{getMethod("enrich.GSE")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/enrich.GSE.R}
#' Users should find it easy to customize this function.
#'
#' @seealso \code{\link{enrich.HGT}}
#' @seealso \code{\link{enrich.DAVID}}
#' @seealso \code{\link{enrich.GOstats}}
#' @seealso \code{\link{enrich.ORT}}
#' @seealso \code{\link{enrichment_analysis}}
#' @seealso \code{\link[DOSE]{enrichResult-class}}
#'
#' @examples
#' data(MLE_Data)
#' universe = TransGeneID(MLE_Data$Gene, "SYMBOL", "ENTREZID", organism = "hsa")
#' geneList = MLE_Data$D7_R1.beta
#' names(geneList) = universe
#' geneList = geneList[!is.na(universe)]
#' enrichRes = enrich.GSE(geneList, type = "KEGG", organism="hsa")
#' head(enrichRes@result)
#'
#'
#' @import clusterProfiler
#' @importFrom pathological temp_dir
#'
#' @export

enrich.GSE <- function(geneList, type= "KEGG", organism='hsa', minGSSize = 10, maxGSSize = 500,
                       pvalueCutoff = 1, pAdjustMethod = "BH"){
  requireNamespace("clusterProfiler", quietly=TRUE) || stop("need clusterProfiler package")
  requireNamespace("pathological", quietly=TRUE) || stop("need pathological package")
  geneList = sort(geneList, decreasing = TRUE)
  #geneList:	order ranked geneList
  if(type == "KEGG"){
    # download Kegg data
    organism = getOrg(organism)$org
    pathwayFiles <- c(file.path(temp_dir(), paste0("pathways_", organism)),
                      file.path(temp_dir(), paste0("gene2path_", organism)))

    if(!all(file.exists(pathwayFiles))){
      gene2path=fread(paste0("http://rest.kegg.jp/link/pathway/",organism),
                      header = FALSE, showProgress = FALSE)
      names(gene2path)=c("EntrezID","PathwayID")
      gene2path$PathwayID=gsub("path:","",gene2path$PathwayID)
      gene2path$EntrezID=gsub(paste0(organism,":"),"",gene2path$EntrezID)

      pathways=fread(paste0("http://rest.kegg.jp/list/pathway/",organism),
                     header = FALSE, showProgress = FALSE)
      names(pathways)=c("PathwayID","PathwayName")

      pathways$PathwayID=gsub("path:","",pathways$PathwayID)
      pathways$PathwayName=gsub(" - .*", "", pathways$PathwayName)

      write.table(pathways, pathwayFiles[1], sep="\t", row.names = FALSE)
      write.table(gene2path, pathwayFiles[2], sep="\t", row.names = FALSE)
    }else{
      pathways=read.table(pathwayFiles[1], sep = "\t", header = TRUE, stringsAsFactors = FALSE)
      gene2path=read.table(pathwayFiles[2], sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    }

    enrichedRes = GSEA(geneList=geneList, minGSSize = minGSSize, maxGSSize = maxGSSize,
                       pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod,
                       TERM2GENE=gene2path[,c("PathwayID","EntrezID")], TERM2NAME=pathways)
  }

  if(type %in% c("BP", "CC", "MF")){
    orgdb = getOrg(organism)$pkg
    enrichedRes = gseGO(geneList=geneList, ont = type, OrgDb=orgdb,
                        minGSSize = minGSSize, maxGSSize = maxGSSize,
                        pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod)
  }

  if(type == "DO"){
    enrichedRes = gseDO(geneList=geneList, minGSSize = minGSSize, maxGSSize = maxGSSize,
                        pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod)
  }
  if(type == "MKEGG"){
    enrichedRes = gseMKEGG(geneList=geneList, organism = organism, minGSSize = minGSSize, maxGSSize = maxGSSize,
                           pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod)
  }
  if(type == "NCG"){
    enrichedRes = gseNCG(geneList=geneList, minGSSize = minGSSize, maxGSSize = maxGSSize,
                           pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod)
  }
  if(!is.null(enrichedRes) && nrow(enrichedRes@result)>0){
    geneID = strsplit(enrichedRes@result$core_enrichment, "/")
    geneName = lapply(geneID, function(gid){
      SYMBOL = TransGeneID(gid, "ENTREZID", "SYMBOL", organism = organism)
      paste(SYMBOL, collapse = "/")
    })
    enrichedRes@result$geneName = unlist(geneName)
  }

  return(enrichedRes)
}

