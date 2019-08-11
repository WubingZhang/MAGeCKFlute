#' Gene set enrichment analysis
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
#' @param type Geneset category for testing, one of 'GOBP', 'GOMF', 'GOCC', 'KEGG',
#' 'BIOCARTA', 'REACTOME', 'CORUM', 'PID', 'HARKMARK', 'c2', 'c6', 'c7', or any
#' combination of them (e.g. 'GOBP+GOMF').
#' @param organism 'hsa' or 'mmu'.
#' @param pvalueCutoff Pvalue cutoff.
#' @param limit A two-length vector (default: c(1, 120)), specifying the minimal and
#' maximal size of gene sets for enrichent analysis.
#' @param gmtpath The path to customized gmt file.
#'
#' @return A enrichResult instance.
#'
#' @author Wubing Zhang
#'
#' @seealso \code{\link{enrich.HGT}}
#' @seealso \code{\link{enrich.ORT}}
#' @seealso \code{\link{EnrichAnalyzer}}
#' @seealso \code{\link[DOSE]{enrichResult-class}}
#'
#' @examples
#' data(geneList, package = "DOSE")
#' \dontrun{
#'     enrichRes = enrich.GSE(geneList)
#'     head(slot(enrichRes, "result"))
#' }
#'
#' @import clusterProfiler
#' @import data.table
#' @import DOSE
#' @export

enrich.GSE <- function(geneList, keytype = "Entrez",
                       type = "CORUM+KEGG",
                       organism = 'hsa', pvalueCutoff = 0.25,
                       limit = c(1, 120), gmtpath = NULL){
  requireNamespace("clusterProfiler", quietly=TRUE) || stop("need clusterProfiler package")
  requireNamespace("data.table", quietly=TRUE) || stop("need data.table package")

  geneList = sort(geneList, decreasing = TRUE)

  ## Prepare gene set annotation
  gene2path = gsGetter(gmtpath, type, limit, organism)
  idx = duplicated(gene2path$PathwayID)
  pathways = data.frame(PathwayID = gene2path$PathwayID[!idx],
                        PathwayName = gene2path$PathwayName[!idx],
                        stringsAsFactors = FALSE)

  ## Gene ID conversion
  if(keytype != "Entrez"){
    allsymbol = names(geneList)
    gene = TransGeneID(allsymbol, keytype, "Entrez", organism = organism)
    idx = is.na(gene) | duplicated(gene)
    geneList = geneList[!idx]
    names(geneList) = gene[!idx]
  }

  ## Enrichment analysis
  len = length(unique(intersect(names(geneList), gene2path$Gene)))
  message("\t", len, " genes are mapped ...")
  enrichedRes = GSEA(geneList = geneList, pvalueCutoff = pvalueCutoff,
                     TERM2NAME = pathways,
                     TERM2GENE = gene2path[,c("PathwayID","Gene")])

  ## Add enriched gene symbols into enrichedRes table
  if(!is.null(enrichedRes) && nrow(enrichedRes@result)>0){
    colnames(enrichedRes@result)[11] = "geneID"
    # enrichedRes@result = EnrichedFilter(enrichedRes@result)
    geneID = strsplit(enrichedRes@result$geneID, "/")
    allsymbol = TransGeneID(names(geneList), "Entrez", "Symbol",
                            organism = organism)
    geneName = lapply(geneID, function(gid){
      SYMBOL = allsymbol[gid]; paste(SYMBOL, collapse = "/")
    })
    enrichedRes@result$geneName = unlist(geneName)
    enrichedRes@result$Count = unlist(lapply(geneID, length))
    cnames = c("ID", "Description", "NES", "pvalue", "p.adjust",
               "geneID", "geneName", "Count")
    enrichedRes@result = enrichedRes@result[, cnames]
  }
  return(enrichedRes)
}

