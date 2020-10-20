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
#' @param type Molecular signatures for testing, available datasets include
#' Pathway (KEGG, REACTOME, C2_CP), GO (GOBP, GOCC, GOMF),
#' MSIGDB (C1, C2 (C2_CP (C2_CP_PID, C2_CP_BIOCARTA), C2_CGP),
#' C3 (C3_MIR, C3_TFT), C4, C6, C7, HALLMARK)
#' and Complex (CORUM). Any combination of them are also accessible
#' (e.g. 'GOBP+GOMF+KEGG+REACTOME').
#' @param method One of "ORT"(Over-Representing Test), "GSEA"(Gene Set Enrichment Analysis), and "HGT"(HyperGemetric test).
#' @param organism 'hsa' or 'mmu'.
#' @param pvalueCutoff FDR cutoff.
#' @param limit A two-length vector (default: c(2, 200)), specifying the minimal and
#' maximal size of gene sets for enrichent analysis.
#' @param universe A character vector, specifying the backgound genelist, default is whole genome.
#' @param filter Boolean, specifying whether filter out redundancies from the enrichment results.
#' @param gmtpath The path to customized gmt file.
#' @param verbose Boolean
#'
#' @return \code{enrichRes} is an enrichResult instance.
#'
#' @author Wubing Zhang
#'
#' @seealso \code{\link{enrich.GSE}}
#' @seealso \code{\link{enrich.ORT}}
#' @seealso \code{\link{enrich.HGT}}
#' @seealso \code{\link[DOSE]{enrichResult-class}}
#'
#' @examples
#' data(geneList, package = "DOSE")
#' \dontrun{
#'   keggA = EnrichAnalyzer(geneList[1:500], keytype = "entrez")
#'   head(keggA@result)
#' }
#' @import stats
#' @export

EnrichAnalyzer = function(geneList, keytype = "Symbol",
                         type = "Pathway+GOBP",
                         method = "HGT",
                         organism = 'hsa',
                         pvalueCutoff = 0.25,
                         limit = c(2, 200),
                         universe = NULL,
                         filter = FALSE,
                         gmtpath = NULL,
                         verbose = TRUE){

  requireNamespace("stats", quietly=TRUE) || stop("need stats package")
  methods = c("ORT", "GSEA", "HGT")
  names(methods) = toupper(methods)
  method = methods[toupper(method)]

  # Gene Set Enrichment Analysis
  if(verbose) message(Sys.time(), " # Running ", type, " enrichment analysis")
  if(method == "GSEA"){
    enrichRes <- enrich.GSE(geneList, keytype = keytype, type = type,
                            organism = organism,
                            pvalueCutoff = pvalueCutoff,
                            limit = limit, gmtpath = gmtpath,
                            verbose = verbose)
  }else if(method == "ORT"){
    enrichRes <- enrich.ORT(geneList, keytype = keytype, type = type,
                            organism = organism, pvalueCutoff = pvalueCutoff,
                            limit = limit, universe = universe,
                            gmtpath = gmtpath, verbose = verbose)
  }else if(method == "HGT"){
    enrichRes = enrich.HGT(geneList, keytype = keytype, type = type,
                           organism = organism, pvalueCutoff = pvalueCutoff,
                           limit = limit, universe = universe,
                           gmtpath = gmtpath, verbose = verbose)
  }else{
    stop("Avaliable methods: GSEA, ORT, and HGT. ")
  }

  if(!is.null(enrichRes) && nrow(enrichRes@result)>10 && filter){
    result = EnrichedFilter(enrichRes)
    # result$p.adjust = p.adjust(result$pvalue, method = "BH")
    # enrichRes@result = result[result$p.adjust<pvalueCutoff, ]
  }
  return(enrichRes)
}

