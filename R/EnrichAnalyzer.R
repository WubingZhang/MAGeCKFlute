#' Enrichment analysis
#'
#' Enrichment analysis
#'
#' @docType methods
#' @name EnrichAnalyzer
#' @rdname EnrichAnalyzer
#' @aliases enrichment
#'
#' @param geneList A numeric vector with gene as names.
#' @param keytype "Entrez" or "Symbol".
#' @param type Geneset category for testing, one of 'CORUM', 'CPX' (ComplexPortal),
#' 'GOBP', 'GOMF', 'GOCC', 'KEGG', 'BIOCARTA', 'REACTOME', 'WikiPathways', 'EHMN', 'PID',
#' or any combination of them (e.g. 'GOBP+GOMF+CORUM'), or 'All' (all categories).
#' @param method One of "ORT"(Over-Representing Test), "GSEA"(Gene Set Enrichment Analysis), and "HGT"(HyperGemetric test).
#' @param organism 'hsa' or 'mmu'.
#' @param pvalueCutoff Pvalue cutoff.
#' @param limit A two-length vector (default: c(3, 50)), specifying the minimal and
#' maximal size of gene sets for enrichent analysis.
#' @param universe A character vector, specifying the backgound genelist, default is whole genome.
#' @param filter Boolean, specifying whether filter out redundancies from the enrichment results.
#' @param gmtpath The path to customized gmt file.
#'
#' @return \code{enrichRes} is an enrichResult instance.
#'
#' @author Feizhen Wu
#'
#' @seealso \code{\link{enrich.GSE}}
#' @seealso \code{\link{enrich.ORT}}
#' @seealso \code{\link{enrich.HGT}}
#' @seealso \code{\link[DOSE]{enrichResult-class}}
#'
#' @examples
#' data(geneList, package = "DOSE")
#' keggA = EnrichAnalyzer(geneList, method = "GSEA")
#' head(keggA@result)
#'
#' @import DOSE
#' @export

EnrichAnalyzer = function(geneList, keytype = "Entrez",
                         type = "CORUM+GOBP+GOMF+GOCC+KEGG",
                         method = "ORT",
                         organism = 'hsa',
                         pvalueCutoff = 0.25,
                         limit = c(3, 80),
                         universe = NULL,
                         filter = TRUE,
                         gmtpath = NA){
  requireNamespace("stats", quietly=TRUE) || stop("need stats package")
  # result = list()
  methods = c("ORT", "GSEA", "HGT")
  names(methods) = toupper(methods)
  method = methods[toupper(method)]

  # Gene Set Enrichment Analysis
  message(Sys.time(), " # Running enrichment analysis for ", type)
  if(method == "GSEA"){
    enrichRes <- enrich.GSE(geneList, keytype = keytype, type = type,
                            organism = organism,
                            pvalueCutoff = 1,
                            limit = limit, gmtpath = gmtpath)
  }else if(method == "ORT"){
    enrichRes <- enrich.ORT(geneList, keytype = keytype, type = type,
                            organism = organism, pvalueCutoff = 1,
                            limit = limit, universe = universe, gmtpath = gmtpath)
  }else if(method == "HGT"){
    enrichRes = enrich.HGT(geneList, keytype = keytype, type = type,
                           organism = organism, pvalueCutoff = 1,
                           limit = limit, universe = universe, gmtpath = gmtpath)
  }else{
    stop("Avaliable methods: GSEA, ORT, and HGT. ")
  }

  if(filter){
    result = EnrichedFilter(enrichRes)
    result$p.adjust = p.adjust(result$pvalue, method = "BH")
    enrichRes@result = result[result$p.adjust<pvalueCutoff, ]
  }
  return(enrichRes)
}

