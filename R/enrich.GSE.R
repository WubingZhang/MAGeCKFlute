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
#' @param type A character, indicating geneset category for testing, "MsigDB_c2_h"(default).
#' @param organism A character, specifying organism, only 'human' is available.
#' @param minGSSize Minimal size of each geneSet for testing.
#' @param maxGSSize Maximal size of each geneSet for analyzing.
#' @param pvalueCutoff Pvalue cutoff.
#' @param pAdjustMethod One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#'
#' @return A enrichResult instance.
#'
#' @author Wubing Zhang
#'
#' @seealso \code{\link{enrich.HGT}}
#' @seealso \code{\link{enrich.DAVID}}
#' @seealso \code{\link{enrich.GOstats}}
#' @seealso \code{\link{enrich.ORT}}
#' @seealso \code{\link{enrichment_analysis}}
#' @seealso \code{\link[DOSE]{enrichResult-class}}
#'
#' @examples
#' data(geneList, package = "DOSE")
#' \dontrun{
#'     enrichRes = enrich.GSE(geneList, type = "KEGG", organism="hsa")
#'     head(enrichRes@result)
#' }
#'
#' @import clusterProfiler
#' @import data.table
#' @import DOSE
#' @export

enrich.GSE <- function(geneList, type = "MsigDB_c2_h", organism='hsa', minGSSize = 10, maxGSSize = 500,
                       pvalueCutoff = 0.25, pAdjustMethod = "BH"){
  requireNamespace("clusterProfiler", quietly=TRUE) || stop("need clusterProfiler package")
  requireNamespace("data.table", quietly=TRUE) || stop("need data.table package")

  geneList = sort(geneList, decreasing = TRUE)
  #geneList:	order ranked geneList
  if(type == "KEGG"){
    # download Kegg data
    organism = getOrg(organism)$org
    pathwayFiles <- c(file.path(system.file("extdata", package = "MAGeCKFlute"),
                                paste0("pathways_", organism)),
                      file.path(system.file("extdata", package = "MAGeCKFlute"),
                                paste0("gene2path_", organism)))
    if(!all(file.exists(pathwayFiles))){
      ## Download pathway annotation
      remfname1 <- paste0("http://rest.kegg.jp/link/pathway/",organism)
      remfname2 <- paste0("http://rest.kegg.jp/list/pathway/",organism)
      download.file(remfname1, pathwayFiles[1], quiet = TRUE)
      download.file(remfname2, pathwayFiles[2], quiet = TRUE)
    }
    ## Read and preprocess pathway annotation
    gene2path = data.table::fread(pathwayFiles[1], header = FALSE, showProgress = FALSE)
    names(gene2path)=c("EntrezID","PathwayID")
    gene2path$PathwayID=gsub("path:","",gene2path$PathwayID)
    gene2path$EntrezID=gsub(paste0(organism,":"),"",gene2path$EntrezID)
    pathways=data.table::fread(pathwayFiles[2], header = FALSE, showProgress = FALSE)
    names(pathways)=c("PathwayID","PathwayName")
    pathways$PathwayID=gsub("path:","",pathways$PathwayID)
    pathways$PathwayName=gsub(" - .*", "", pathways$PathwayName)
    ##==========
    enrichedRes = GSEA(geneList=geneList, minGSSize = minGSSize, maxGSSize = maxGSSize,
                       pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod,
                       TERM2GENE=gene2path[,c("PathwayID","EntrezID")], TERM2NAME=pathways)
  }

  if(type == "MsigDB_c2_h"){
    gene2path = read.gmt(system.file("extdata", "MsigDB_c2_h.gmt", package = "MAGeCKFlute"))
    enrichedRes = GSEA(geneList=geneList, minGSSize = minGSSize, maxGSSize = maxGSSize,
                       pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod,
                       TERM2GENE=gene2path)
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
    allsymbol = TransGeneID(names(geneList), "Entrez", "Symbol", organism = organism)
    geneName = lapply(geneID, function(gid){
      SYMBOL = allsymbol[gid]
      paste(SYMBOL, collapse = "/")
    })
    enrichedRes@result$geneName = unlist(geneName)
  }

  return(enrichedRes)
}

