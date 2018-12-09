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
#' @param universe A character vector, specifying the backgound genelist, default is whole genome.
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
#' @return A enrichedResult instance.
#'
#' @author Wubing Zhang
#'
#' @seealso \code{\link{enrich.HGT}}
#' @seealso \code{\link{enrich.GSE}}
#' @seealso \code{\link{enrichment_analysis}}
#' @seealso \code{\link[DOSE]{enrichResult-class}}
#'
#' @examples
#' data(geneList, package = "DOSE")
#' genes <- geneList[1:100]
#' enrichedRes <- enrich.ORT(genes)
#' head(enrichedRes@result)
#'
#' @import DOSE
#' @import clusterProfiler
#' @export

enrich.ORT <- function(geneList, universe=NULL, keytype = "Entrez", type="GOBP+GOMF", organism='hsa',
                       pvalueCutoff = 0.25, pAdjustMethod = "BH", limit = c(3, 50), gmtpath = NA){
  requireNamespace("clusterProfiler", quietly=TRUE) || stop("need clusterProfiler package")
  requireNamespace("data.table", quietly=TRUE) || stop("need data.table package")

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
                        PathwayName = gene2path$PathwayName[!idx],
                        stringsAsFactors = FALSE)

  ## Gene ID conversion
  if(keytype == "Symbol"){
    allsymbol = names(geneList)
    gene = TransGeneID(allsymbol, "Symbol", "Entrez", organism = organism)
    geneList = geneList[!duplicated(gene)]
    names(geneList) = gene[!duplicated(gene)]
    universe = TransGeneID(universe, "Symbol", "Entrez", organism = organism)
    universe = universe[!is.na(universe)]
  }
  gene = names(geneList)
  gene2path$Gene = as.character(gene2path$Gene)
  if(!is.null(universe)){
    universe = unique(as.character(universe))
  }else{
    universe = unique(c(gene, gene2path$Gene))
  }

  ## Enrichment analysis
  orgdb = getOrg(organism)$pkg
  enrichedRes = enricher(gene, pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod,
                universe = universe, minGSSize = limit[1], maxGSSize = limit[2],
                TERM2GENE = gene2path[,c("PathwayID","Gene")], TERM2NAME = pathways)

  ## Add enriched gene symbols into enrichedRes table
  if(!is.null(enrichedRes)){
    allsymbol = TransGeneID(gene, "Entrez", "Symbol", organism = organism)
    geneID = strsplit(enrichedRes@result$geneID, split = "/")
    geneName = lapply(geneID, function(gid){
      SYMBOL = allsymbol[gid]
      paste(SYMBOL, collapse = "/")
    })
    enrichedRes@result$geneName = unlist(geneName)
    enrichedRes@result$NES = unlist(apply(enrichedRes@result, 1, function(x){
      enrichGenes = unlist(strsplit(x["geneID"], split = "/"))
      NES = sum(geneList[enrichGenes]) / log2(length(enrichGenes)+1)
      return(NES)
    }))
    idx = c("ID", "Description", "NES", "pvalue", "p.adjust",
            "GeneRatio", "BgRatio", "geneID", "geneName", "Count")
    enrichedRes@result = enrichedRes@result[, idx]
  }
  return(enrichedRes)
}

