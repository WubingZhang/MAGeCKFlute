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
#' @return A list containing enrichment results for each group genes. This list contains eight
#' items, which contain subitems of \code{gridPlot} and \code{enrichRes}.
#'
#' @author Wubing Zhang
#'
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
  idx1 = gg$group=="top"; genes = gg$EntrezID[idx1]
  geneList = gg$Diff[idx1]; names(geneList) = genes
  enrichA = EnrichAnalyzer(geneList = geneList, universe = gg$EntrezID,
                         method = enrich_method, type = "KEGG+REACTOME+GOBP+Complex",
                         organism = organism, pvalueCutoff = pvalue,
                         limit = limit, keytype = "Entrez")
  if(!is.null(enrichA) && nrow(enrichA@result)>0){
    keggA = enrichA@result[grepl("KEGG", enrichA@result$ID), ]
    gobpA = enrichA@result[grepl("GOBP", enrichA@result$ID), ]
    reactomeA = enrichA@result[grepl("REACTOME", enrichA@result$ID), ]
    complexA = enrichA@result[grepl("CPX|CORUM", enrichA@result$ID), ]
    keggA = list(enrichRes = keggA, gridPlot = EnrichedView(keggA, top = 5, bottom = 0)
                 + labs(title = "KEGG: GroupA"))
    gobpA = list(enrichRes = gobpA, gridPlot = EnrichedView(gobpA, top = 5, bottom = 0)
               + labs(title = "GOBP: GroupA"))
    reactomeA = list(enrichRes = reactomeA, gridPlot = EnrichedView(reactomeA, top = 5, bottom = 0)
                     + labs(title = "REACTOME: GroupA"))
    complexA = list(enrichRes = complexA, gridPlot = EnrichedView(complexA, top = 5, bottom = 0)
                    + labs(title = "Complex: GroupA"))
  }else{
    keggA = gobpA = reactomeA = complexA = list(enrichRes = NULL, gridPlot = noEnrichPlot())
  }

  idx2 = gg$group=="bottom"; genes = gg$EntrezID[idx2]
  geneList = gg$Diff[idx2]; names(geneList) = genes
  enrichB = EnrichAnalyzer(geneList = geneList, universe = gg$EntrezID,
                           method = enrich_method, type = "KEGG+REACTOME+GOBP+Complex",
                           organism = organism, pvalueCutoff = pvalue,
                           limit = limit, keytype = "Entrez")
  if(!is.null(enrichB) && nrow(enrichB@result)>0){
    keggB = enrichB@result[grepl("KEGG", enrichB@result$ID), ]
    gobpB = enrichB@result[grepl("GOBP", enrichB@result$ID), ]
    reactomeB = enrichB@result[grepl("REACTOME", enrichB@result$ID), ]
    complexB = enrichB@result[grepl("CPX|CORUM", enrichB@result$ID), ]
    keggB = list(enrichRes = keggB, gridPlot = EnrichedView(keggB, top = 0, bottom = 5)
                 + labs(title = "KEGG: GroupB"))
    gobpB = list(enrichRes = gobpB, gridPlot = EnrichedView(gobpB, top = 0, bottom = 5)
               + labs(title = "GOBP: GroupB"))
    reactomeB = list(enrichRes = reactomeB, gridPlot = EnrichedView(reactomeB, top = 0, bottom = 5)
                     + labs(title = "REACTOME: GroupB"))
    complexB = list(enrichRes = complexB, gridPlot = EnrichedView(complexB, top = 0, bottom = 5)
                    + labs(title = "Complex: GroupB"))
  }else{
    keggB = gobpB = reactomeB = complexB = list(enrichRes = NULL, gridPlot = noEnrichPlot())
  }

  if(!is.null(filename)){
    ## Save GroupA enrichment results
    if(!is.null(enrichA) && nrow(enrichA@result)>0){
      write.table(keggA$enrichRes, file.path(out.dir,paste0("GroupA_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(keggA$gridPlot, filename=file.path(out.dir,paste0("GroupA_kegg_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(reactomeA$enrichRes, file.path(out.dir,paste0("GroupA_reactome_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(reactomeA$gridPlot, filename=file.path(out.dir,paste0("GroupA_reactome_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(gobpA$enrichRes, file.path(out.dir,paste0("GroupA_gobp_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(gobpA$gridPlot, filename=file.path(out.dir,paste0("GroupA_gobp_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(complexA$enrichRes, file.path(out.dir,paste0("GroupA_complex_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(complexA$gridPlot, filename=file.path(out.dir,paste0("GroupA_complex_", filename,".png")),
             units = "in", width=6.5, height=4)
    }
    if(!is.null(enrichB) && nrow(enrichB@result)>0){
      write.table(keggB$enrichRes, file.path(out.dir,paste0("GroupB_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(keggB$gridPlot, filename=file.path(out.dir,paste0("GroupB_kegg_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(reactomeB$enrichRes, file.path(out.dir,paste0("GroupB_reactome_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(reactomeB$gridPlot, filename=file.path(out.dir,paste0("GroupB_reactome_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(gobpB$enrichRes, file.path(out.dir,paste0("GroupB_gobp_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(gobpB$gridPlot, filename=file.path(out.dir,paste0("GroupB_gobp_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(complexB$enrichRes, file.path(out.dir,paste0("GroupB_complex_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(complexB$gridPlot, filename=file.path(out.dir,paste0("GroupB_complex_", filename,".png")),
             units = "in", width=6.5, height=4)
    }
  }
  return(list(keggA=keggA, gobpA=gobpA, reactomeA=reactomeA, complexA=complexA,
              keggB=keggB, gobpB=gobpB, reactomeB=reactomeB, complexB=complexB))
}
