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
#' @param type Molecular signatures for testing, available datasets include
#' Pathway (PID, KEGG, REACTOME, BIOCARTA, C2CP), GO (GOBP, GOCC, GOMF),
#' Complex (CORUM, CPX), c1, c2, c3, c4, c6, c7, HALLMARK. It also allows any
#' combination of them (e.g. 'GOBP+GOMF+KEGG+REACTOME').
#' @param organism 'hsa' or 'mmu'.
#' @param pvalueCutoff Pvalue cutoff.
#' @param limit A two-length vector, specifying the minimal and
#' maximal size of gene sets for enrichent analysis.
#' @param gmtpath The path to customized gmt file.
#' @param nPerm The number of permutations.
#' @param by One of 'fgsea' or 'DOSE'
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
#'     enrichRes = enrich.GSE(geneList, keytype = "entrez")
#'     head(slot(enrichRes, "result"))
#' }
#'
#' @import data.table DOSE
#' @export

enrich.GSE <- function(geneList, keytype = "Symbol",
                       type = "Pathway+GOBP",
                       organism = 'hsa', pvalueCutoff = 0.25,
                       limit = c(2, 200), gmtpath = NULL,
                       nPerm = 2000, by = "fgsea"){
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
                     minGSSize = 0, maxGSSize = max(limit),
                     TERM2NAME = pathways, nPerm = nPerm, by = by,
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
    res = enrichedRes@result[, cnames]
    res = res[order(res$pvalue), ]
    enrichedRes@result = res
  }
  return(enrichedRes)
}

