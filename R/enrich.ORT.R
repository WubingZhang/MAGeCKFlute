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
#' @param limit A two-length vector (default: c(3, 50)), specifying the minimal and
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
                       type = "CORUM+GOBP+GOMF+GOCC+KEGG",
                       organism = 'hsa', pvalueCutoff = 0.05,
                       limit = c(3, 80), universe=NULL, gmtpath = NA){
  requireNamespace("clusterProfiler", quietly=TRUE) || stop("need clusterProfiler package")
  requireNamespace("data.table", quietly=TRUE) || stop("need data.table package")

  ## Prepare gene set annotation
  if(is.na(gmtpath)){
    msigdb = file.path(system.file("extdata", package = "MAGeCKFlute"),
                       paste0(organism, "_msig_entrez.gmt.gz"))
    gmtpath = gzfile(msigdb)
  }
  gene2path = ReadGMT(gmtpath, limit = limit)
  close(gmtpath)
  names(gene2path) = c("Gene","PathwayID", "PathwayName")
  gene2path$PathwayName = paste0(toupper(substr(gene2path$PathwayName, 0, 1)),
                                 tolower(substr(gene2path$PathwayName, 2,
                                                nchar(gene2path$PathwayName))))

  if(type != "All"){
    type = unlist(strsplit(type, "\\+"))
    idx = toupper(gsub("_.*", "", gene2path$PathwayID)) %in% toupper(type)
    gene2path = gene2path[idx, ]
  }
  gene2path = gene2path[!is.na(gene2path$Gene), ]
  idx = duplicated(gene2path$PathwayID)
  pathways = data.frame(PathwayID = gene2path$PathwayID[!idx],
                        PathwayName = gene2path$PathwayName[!idx],
                        stringsAsFactors = FALSE)

  ## Gene ID conversion
  if(keytype != "Entrez"){
    allsymbol = names(geneList)
    gene = TransGeneID(allsymbol, keytype, "Entrez", organism = organism)
    geneList = geneList[!duplicated(gene)]
    names(geneList) = gene[!duplicated(gene)]
    universe = TransGeneID(universe, keytype, "Entrez", organism = organism)
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

