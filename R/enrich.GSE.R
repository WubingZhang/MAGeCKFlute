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
#' @param keytype "Entrez" or "Symbol".
#' @param type Geneset category for testing, "MsigDB" (default), "KEGG" or "Others" (`gmtpath` is required).
#' @param organism A character, specifying organism, only 'human' is available.
#' @param pvalueCutoff Pvalue cutoff.
#' @param pAdjustMethod One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param limit A integer vector of length two, specifying the limit of geneset size.
#' @param gmtpath The path to gmt file.
#'
#' @return A enrichResult instance.
#'
#' @author Wubing Zhang
#'
#' @seealso \code{\link{enrich.HGT}}
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

enrich.GSE <- function(geneList, keytype = "Entrez", type="MsigDB", organism='hsa',
                       pvalueCutoff = 0.05, pAdjustMethod = "BH", limit = c(3, 50), gmtpath = NA){
  requireNamespace("clusterProfiler", quietly=TRUE) || stop("need clusterProfiler package")
  requireNamespace("data.table", quietly=TRUE) || stop("need data.table package")

  geneList = sort(geneList, decreasing = TRUE)
  ## Gene ID conversion
  if(keytype == "Symbol"){
    allsymbol = names(geneList)
    gene = TransGeneID(allsymbol, "Symbol", "Entrez", organism = organism)
    geneList = geneList[!duplicated(gene)]
    names(geneList) = gene[!duplicated(gene)]
  }

  ## Enrichment analysis on MsigDB genesets
  if(type %in% c("MsigDB", "Others")){
    msigdb = file.path(system.file("extdata", package = "MAGeCKFlute"),
                       paste0("MsigDB.all.v6.2.entrez.", organism, ".gmt.gz"))
    if(type == "MsigDB"){ gmtpath = gzfile(msigdb) }
    gene2path = readGMT(gmtpath, limit = limit)
    names(gene2path) = c("Gene","PathwayID", "PathwayName")
    gene2path = gene2path[!is.na(gene2path$Gene), ]
    pathways = data.frame(PathwayID = unique(gene2path$PathwayID),
                          PathwayName = unique(gene2path$PathwayName))
    enrichedRes = GSEA(geneList = geneList, minGSSize = limit[1], maxGSSize = limit[2],
                       pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod,
                       TERM2GENE = gene2path[,c("PathwayID","Gene")], TERM2NAME = pathways)
  }

  if(type == "KEGG"){
    # download Kegg data
    organism = getOrg(organism)$org
    pathwayFiles <- c(file.path(system.file("extdata", package = "MAGeCKFlute"),
                                paste0("pathways_", organism)),
                      file.path(system.file("extdata", package = "MAGeCKFlute"),
                                paste0("gene2path_", organism)))
    if(!all(file.exists(pathwayFiles))){
      ## Download pathway annotation
      remfname1 <- paste0("http://rest.kegg.jp/link/pathway/", organism)
      remfname2 <- paste0("http://rest.kegg.jp/list/pathway/", organism)
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
    enrichedRes = GSEA(geneList=geneList, minGSSize = limit[1], maxGSSize = limit[2],
                       pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod,
                       TERM2GENE=gene2path[,c("PathwayID","EntrezID")], TERM2NAME=pathways)
  }

  if(type %in% c("BP", "CC", "MF")){
    orgdb = getOrg(organism)$pkg
    enrichedRes = gseGO(geneList=geneList, ont = type, OrgDb=orgdb,
                        minGSSize = limit[1], maxGSSize = limit[2],
                        pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod)
  }

  if(type == "DO"){
    enrichedRes = gseDO(geneList=geneList, minGSSize = limit[1], maxGSSize = limit[2],
                        pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod)
  }

  if(type == "MKEGG"){
    enrichedRes = gseMKEGG(geneList=geneList, organism = organism, minGSSize = limit[1], maxGSSize = limit[2],
                           pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod)
  }

  if(type == "NCG"){
    enrichedRes = gseNCG(geneList=geneList, minGSSize = limit[1], maxGSSize = limit[2],
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

