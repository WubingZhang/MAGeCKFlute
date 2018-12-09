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
#' @param type Geneset category for testing, one of 'GOBP+GOMF' (default), 'GOBP', 'GOMF', 'GOCC',
#' 'KEGG', 'BIOCARTA', 'REACTOME', 'WikiPathways', 'EHMN', 'PID', or 'All' and any combination of them,
#' such as 'KEGG+BIOCARTA+REACTOME+GOBP+GOCC+GOMF+EHMN+PID+WikiPathways'.
#' @param organism 'hsa' or 'mmu'.
#' @param pvalueCutoff Pvalue cutoff.
#' @param pAdjustMethod One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param limit A two-length vector (default: c(3, 50)), specifying the minimal and
#' maximal size of gene sets for enrichent analysis.
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

enrich.GSE <- function(geneList, keytype = "Entrez", type="GOBP+GOMF", organism='hsa',
                       pvalueCutoff = 0.25, pAdjustMethod = "BH", limit = c(3, 50), gmtpath = NA){
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

  ## Prepare gene set annotation
  if(is.na(gmtpath)){
    msigdb = file.path(system.file("extdata", package = "MAGeCKFlute"),
                       paste0(organism, "_msig_entrez.gmt.gz"))
    gmtpath = gzfile(msigdb)
  }
  gene2path = ReadGMT(gmtpath, limit = c(1, 500))
  close(gmtpath)
  names(gene2path) = c("Gene","PathwayID", "PathwayName")
  gene2path$PathwayName = toupper(gene2path$PathwayName)
  if(type == "All") type = 'KEGG+BIOCARTA+REACTOME+GOBP+GOCC+GOMF+EHMN+PID+WikiPathways'
  type = unlist(strsplit(type, "\\+"))
  idx = toupper(gsub("_.*", "", gene2path$PathwayID)) %in% toupper(type)
  gene2path = gene2path[idx, ]
  gene2path = gene2path[!is.na(gene2path$Gene), ]
  idx = duplicated(gene2path$PathwayID)
  pathways = data.frame(PathwayID = gene2path$PathwayID[!idx],
                        PathwayName = gene2path$PathwayName[!idx])

  ## Enrichment analysis
  enrichedRes = GSEA(geneList = geneList, minGSSize = limit[1], maxGSSize = limit[2],
                     pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod,
                     TERM2GENE = gene2path[,c("PathwayID","Gene")], TERM2NAME = pathways)

  ## Add enriched gene symbols into enrichedRes table
  if(!is.null(enrichedRes) && nrow(enrichedRes@result)>0){
    geneID = strsplit(enrichedRes@result$core_enrichment, "/")
    allsymbol = TransGeneID(names(geneList), "Entrez", "Symbol", organism = organism)
    geneName = lapply(geneID, function(gid){
      SYMBOL = allsymbol[gid]
      paste(SYMBOL, collapse = "/")
    })
    enrichedRes@result$geneName = unlist(geneName)
    enrichedRes@result$Count = unlist(lapply(geneID, length))
    enrichedRes@result = enrichedRes@result[, c("ID", "Description", "NES", "pvalue", "p.adjust",
                                                "core_enrichment", "geneName", "Count")]
    colnames(enrichedRes@result)[6] = "geneID"
  }

  return(enrichedRes)
}

