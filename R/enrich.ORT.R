#' Do enrichment analysis using over-representation test
#'
#' Do enrichment analysis using over-representation test
#'
#' @docType methods
#' @name enrich.ORT
#' @rdname enrich.ORT
#' @aliases enrichORT
#'
#' @param geneList A numeric vector with gene as names.
#' @param keytype "Entrez" or "Symbol".
#' @param type Geneset category for testing, one of 'CORUM', 'CPX' (ComplexPortal),
#' 'GOBP', 'GOMF', 'GOCC', 'KEGG', 'BIOCARTA', 'REACTOME', 'WikiPathways', 'EHMN', 'PID',
#' or any combination of them (e.g. 'GOBP+GOMF+CORUM'), or 'All' (all categories).
#' @param organism 'hsa' or 'mmu'.
#' @param pvalueCutoff Pvalue cutoff.
#' @param limit A two-length vector (default: c(1, 120)), specifying the minimal and
#' maximal size of gene sets for enrichent analysis.
#' @param universe A character vector, specifying the backgound genelist, default is whole genome.
#' @param gmtpath The path to customized gmt file.
#'
#' @return A enrichedResult instance.
#'
#' @author Wubing Zhang
#'
#' @seealso \code{\link{enrich.HGT}}
#' @seealso \code{\link{enrich.GSE}}
#' @seealso \code{\link{EnrichAnalyzer}}
#' @seealso \code{\link[DOSE]{enrichResult-class}}
#'
#' @examples
#' data(geneList, package = "DOSE")
#' genes <- geneList[1:100]
#' enrichedRes <- enrich.ORT(genes)
#' head(slot(enrichedRes, "result"))
#'
#' @import DOSE
#' @import clusterProfiler
#' @export

enrich.ORT <- function(geneList, keytype = "Entrez",
                       type = "CORUM+KEGG",
                       organism = 'hsa', pvalueCutoff = 0.05,
                       limit = c(1, 120), universe=NULL, gmtpath = NA){
  requireNamespace("clusterProfiler", quietly=TRUE) || stop("need clusterProfiler package")
  requireNamespace("data.table", quietly=TRUE) || stop("need data.table package")

  ## Prepare gene set annotation
  gene2path = gsGetter(gmtpath, type, limit, organism)
  idx = duplicated(gene2path$PathwayID)
  pathways = data.frame(PathwayID = gene2path$PathwayID[!idx],
                        PathwayName = gene2path$PathwayName[!idx],
                        stringsAsFactors = FALSE)

  ## Gene ID conversion
  gene = names(geneList)
  if(keytype != "Entrez"){
    allsymbol = names(geneList)
    gene = TransGeneID(allsymbol, keytype, "Entrez", organism = organism)
    idx = duplicated(gene)|is.na(gene)
    geneList = geneList[!idx]; names(geneList) = gene[!idx]
  }
  if(!is.null(universe)){
    universe = TransGeneID(universe, keytype, "Entrez", organism = organism)
    universe = universe[!is.na(universe)]
  }else{
    universe = unique(c(gene, gene2path$Gene))
  }
  gene = names(geneList)

  ## Enrichment analysis
  len = length(unique(intersect(gene, gene2path$Gene)))
  message("\t", len, " genes are mapped ...")
  orgdb = getOrg(organism)$pkg
  enrichedRes = enricher(gene, universe = universe,
                         TERM2NAME = pathways, pvalueCutoff = pvalueCutoff,
                         TERM2GENE = gene2path[,c("PathwayID","Gene")])
  if(!is.null(enrichedRes) && nrow(enrichedRes@result)>0){
    enrichedRes@result = enrichedRes@result[enrichedRes@result$p.adjust<pvalueCutoff, ]
  }
  ## Add enriched gene symbols into enrichedRes table
  if(!is.null(enrichedRes) && nrow(enrichedRes@result)>0){
    allsymbol = TransGeneID(gene, "Entrez", "Symbol", organism = organism)
    geneID = strsplit(enrichedRes@result$geneID, split = "/")
    geneName = lapply(geneID, function(gid){
      SYMBOL = allsymbol[gid]; paste(SYMBOL, collapse = "/")
    })
    enrichedRes@result$geneName = unlist(geneName)
    enrichedRes@result$NES = as.vector(sapply(enrichedRes@result$geneID, function(x){
      enrichGenes = unlist(strsplit(x, split = "/"))
      NES = mean(geneList[enrichGenes]) * length(enrichGenes)^0.6
      return(NES)
    }))
    idx = c("ID", "Description", "NES", "pvalue", "p.adjust",
            "GeneRatio", "BgRatio", "geneID", "geneName", "Count")
    enrichedRes@result = enrichedRes@result[, idx]
  }
  return(enrichedRes)
}

