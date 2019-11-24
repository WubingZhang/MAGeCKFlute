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
#' @param type Molecular signatures for testing, available datasets include
#' Pathway (PID, KEGG, REACTOME, BIOCARTA, C2CP), GO (GOBP, GOCC, GOMF),
#' Complex (CORUM, CPX), c1, c2, c3, c4, c6, c7, HALLMARK. It also allows any
#' combination of them (e.g. 'GOBP+GOMF+KEGG+REACTOME').
#' @param organism 'hsa' or 'mmu'.
#' @param pvalueCutoff Pvalue cutoff.
#' @param limit A two-length vector, specifying the minimal and
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
#' enrichedRes <- enrich.ORT(genes, keytype = "entrez")
#' head(slot(enrichedRes, "result"))
#'
#' @import DOSE
#' @export

enrich.ORT <- function(geneList, keytype = "Symbol",
                       type = "Pathway+GOBP",
                       organism = 'hsa', pvalueCutoff = 0.25,
                       limit = c(2, 200), universe=NULL, gmtpath = NULL){
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
                         minGSSize = 0, maxGSSize = max(limit),
                         TERM2NAME = pathways, pvalueCutoff = pvalueCutoff,
                         TERM2GENE = gene2path[,c("PathwayID","Gene")])
  if(!is.null(enrichedRes) && nrow(enrichedRes@result)>0){
    res = enrichedRes@result[enrichedRes@result$p.adjust<=pvalueCutoff, ]
    res = res[order(res$pvalue), ]
    enrichedRes@result = res
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

