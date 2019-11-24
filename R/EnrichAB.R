#' Enrichment analysis for Positive and Negative selection genes
#'
#' Do enrichment analysis for selected genes, in which positive selection and negative selection
#' are termed as GroupA and GroupB
#'
#' @docType methods
#' @name EnrichAB
#' @rdname EnrichAB
#'
#' @param data A data frame.
#' @param pvalue Pvalue cutoff.
#' @param enrich_method One of "ORT"(Over-Representing Test), "GSEA"(Gene Set Enrichment Analysis), and "HGT"(HyperGemetric test).
#' @param organism "hsa" or "mmu".
#' @param limit A two-length vector (default: c(1, 120)), specifying the min and
#' max size of pathways for enrichent analysis.
#' @param filename Suffix of output file name.
#' @param out.dir Path to save plot to (combined with filename).
#' @param width As in ggsave.
#' @param height As in ggsave.
#' @param ... Other available parameters in ggsave.
#'
#' @return A list containing enrichment results for each group genes. This list contains items four
#' items, \code{keggA}, \code{keggB}, \code{goA}, \code{goB}. Four items are all list object, containing
#' subitems of \code{gridPlot} and \code{enrichRes}. \code{gridPlot} is a ggplot object, and
#' \code{enrichRes} is a enrichResult instance
#'
#' @author Wubing Zhang

# Enrichment for GroupA and GroupB genes
EnrichAB <- function(data, pvalue = 0.25,
                     enrich_method = "ORT",
                     organism = "hsa",
                     limit = c(1, 120),
                     filename = NULL, out.dir = ".",
                     width = 6.5, height = 4, ...){

  requireNamespace("clusterProfiler", quietly=TRUE) || stop("Need clusterProfiler package")
  message(Sys.time(), " # Enrichment analysis of GroupA and GroupB genes ...")
  gg = data

  ##=====enrichment for GroupA======
  idx1 = gg$group=="up"; genes = gg$EntrezID[idx1]
  geneList = gg$diff[idx1]; names(geneList) = genes
  keggA = EnrichAnalyzer(geneList = geneList, universe = gg$EntrezID,
                         method = enrich_method, type = "KEGG",
                         organism = organism, pvalueCutoff = pvalue,
                         limit = limit)
  keggA = list(enrichRes = keggA,
               gridPlot = EnrichedView(keggA, top = 8, bottom = 0)
               + labs(title = "KEGG: GroupA"))
  goA = EnrichAnalyzer(geneList = geneList, universe = gg$EntrezID,
                       method = "ORT", type = "GOBP+GOMF+GOCC",
                       organism = organism, pvalueCutoff = pvalue,
                       limit = limit)
  goA = list(enrichRes = goA,
             gridPlot = EnrichedView(goA, top = 8, bottom = 0)
               + labs(title = "Gene Ontology: GroupA"))


  idx2 = gg$group=="down"; genes = gg$EntrezID[idx2]
  geneList = gg$diff[idx2]; names(geneList) = genes
  keggB=EnrichAnalyzer(geneList = geneList, universe = gg$EntrezID,
                       method = enrich_method, type = "KEGG",
                       organism = organism, pvalueCutoff = pvalue,
                       limit = limit)
  keggB = list(enrichRes = keggB,
               gridPlot = EnrichedView(keggB, top = 0, bottom = 8)
               + labs(title = "KEGG: GroupB"))
  goB = EnrichAnalyzer(geneList = geneList, universe = gg$EntrezID,
                       method = "ORT", type = "GOBP+GOMF+GOCC",
                       organism = organism, pvalueCutoff = pvalue,
                       limit = limit)
  goB = list(enrichRes = goB,
             gridPlot = EnrichedView(goB, top = 0, bottom = 8)
               + labs(title = "Gene Ontology: GroupB"))

  if(!is.null(filename)){
    ## Save GroupA enrichment results
    if(!is.null(keggA$enrichRes)){
      write.table(keggA$enrichRes@result,
                  file.path(out.dir,paste0("GroupA_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(keggA$gridPlot,
             filename=file.path(out.dir,paste0("GroupA_kegg_", filename,".png")),
             units = "in", width=6.5, height=4)
    }
    if(!is.null(goA$enrichRes)){
      write.table(goA$enrichRes@result,
                  file.path(out.dir,paste0("GroupA_go_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(goA$gridPlot,
             filename=file.path(out.dir,paste0("GroupA_go_",filename,".png")),
             units = "in", width=6.5, height=4)
    }
    ##=========Save GroupB enrichment results===========================
    if(!is.null(keggB$enrichRes)){
      write.table(keggB$enrichRes@result,
                  file.path(out.dir,paste0("GroupB_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(keggB$gridPlot,
             filename=file.path(out.dir,paste0("GroupB_kegg_",filename,".png")),
             units = "in", width=6.5, height=4)
    }
    if(!is.null(goB$enrichRes)){
      write.table(goB$enrichRes@result,
                  file.path(out.dir,paste0("GroupB_go_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
      ggsave(goB$gridPlot,
             filename=file.path(out.dir,paste0("GroupB_go_",filename,".png")),
             units = "in", width=6.5, height=4)
    }
    return(list(keggA=keggA, goA=goA, keggB=keggB, goB=goB))
  }
}
