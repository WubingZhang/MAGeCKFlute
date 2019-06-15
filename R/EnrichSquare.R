#' Enrichment analysis for selected treatment related genes
#'
#' Do enrichment analysis for selected treatment related genes in 9-squares
#'
#' @docType methods
#' @name EnrichSquare
#' @rdname EnrichSquare
#'
#' @param beta Data frame, with rownames of Entrez IDs, which contains columns of 'group' and 'diff'.
#' @param pvalue Pvalue cutoff.
#' @param enrich_method One of "ORT"(Over-Representing Test) and "HGT"(HyperGemetric test).
#' @param organism "hsa" or "mmu".
#' @param limit A two-length vector (default: c(3, 50)), specifying the min and
#' max size of pathways for enrichent analysis.
#' @param filename Suffix of output file name. NULL(default) means no output.
#' @param out.dir Path to save plot to (combined with filename).
#' @param width As in ggsave.
#' @param height As in ggsave.
#' @param ... Other available parameters in ggsave.
#'
#' @return A list containing enrichment results for each group genes. This list contains several elements:
#' \item{kegg1}{a list record enriched KEGG pathways for Group1 genes in 9-Square}
#' \item{kegg2}{a list record enriched KEGG pathways for Group2 genes in 9-Square}
#' \item{kegg3}{a list record enriched KEGG pathways for Group3 genes in 9-Square}
#' \item{kegg4}{a list record enriched KEGG pathways for Group4 genes in 9-Square}
#' \item{kegg12}{a list record enriched KEGG pathways for Group1&Group2 genes in 9-Square}
#' \item{kegg13}{a list record enriched KEGG pathways for Group1&Group3 genes in 9-Square}
#' \item{kegg24}{a list record enriched KEGG pathways for Group2&Group4 genes in 9-Square}
#' \item{kegg34}{a list record enriched KEGG pathways for Group3&Group4 genes in 9-Square}
#' \item{go1}{a list record enriched GO terms for Group1 genes in 9-Square}
#' \item{go2}{a list record enriched GO terms for Group2 genes in 9-Square}
#' \item{go3}{a list record enriched GO terms for Group3 genes in 9-Square}
#' \item{go4}{a list record enriched GO terms for Group4 genes in 9-Square}
#' \item{go12}{a list record enriched GO terms for Group1&Group2 genes in 9-Square}
#' \item{go13}{a list record enriched GO terms for Group1&Group3 genes in 9-Square}
#' \item{go24}{a list record enriched GO terms for Group2&Group4 genes in 9-Square}
#' \item{go34}{a list record enriched GO terms for Group3&Group4 genes in 9-Square}
#'
#' Each item in the returned list has two sub items:
#' \item{gridPlot}{an object created by \code{ggplot}, which can be assigned and further customized.}
#' \item{enrichRes}{a enrichResult instance.}
#'
#' @author Wubing Zhang
#'
#'
#' @seealso \code{\link{SquareView}}
#' @seealso \code{\link{EnrichSquare}}
#'
#' @examples
#' data(mle.gene_summary)
#' dd = ReadBeta(mle.gene_summary)
#' p = SquareView(dd, ctrlname = "dmso", treatname = "plx")
#' \dontrun{
#'  # Read beta score from gene summary table in MAGeCK MLE results
#'  E1 = EnrichSquare(p$data, organism="hsa")
#'  print(E1$kegg1$gridPlot)
#'}
#' @export

#enrichment for square grouped genes
EnrichSquare <- function(beta, pvalue = 0.05, enrich_method="ORT", organism="hsa", limit = c(3,50),
                         filename=NULL, out.dir=".", width=6.5, height=4, ...){
  message(Sys.time(), " # Enrichment analysis of 9 Square grouped genes ...")

  ## ===========Enrichment===================
  gg=beta
  idx1 = gg$group=="midleft"
  idx2 = gg$group=="topcenter"
  idx3 = gg$group=="midright"
  idx4 = gg$group=="bottomcenter"
  universe = rownames(gg)

  #====GO_KEGG_enrichment=====
  geneList = gg$diff[idx1]; names(geneList) = rownames(gg)[idx1]
  kegg1=EnrichAnalyzer(geneList = geneList, universe=universe,
                       method = enrich_method, type = "KEGG",
                       organism=organism, pvalueCutoff = pvalue,
                       limit = limit)
  go1=EnrichAnalyzer(geneList = geneList, universe=universe,
                     method = enrich_method, type = "GOBP+GOMF+GOCC",
                     organism=organism, pvalueCutoff = pvalue,
                     limit = limit)
  kegg1 = list(enrichRes = kegg1,
               gridPlot = EnrichedView(EnrichedFilter(kegg1), top = 8, bottom = 0)
               + labs(title = "KEGG: Group1"))
  go1 = list(enrichRes = go1,
             gridPlot = EnrichedView(EnrichedFilter(go1), top = 8, bottom = 0)
               + labs(title = "Gene Ontology: Group1"))


  #====GO_KEGG_enrichment=====
  geneList = gg$diff[idx2]; names(geneList) = rownames(gg)[idx2]
  kegg2=EnrichAnalyzer(geneList = geneList, universe=universe,
                       method = enrich_method, type = "KEGG",
                       organism=organism, pvalueCutoff=pvalue,
                       limit = limit)
  go2=EnrichAnalyzer(geneList = geneList, universe=universe,
                     method = enrich_method, type = "GOBP+GOMF+GOCC",
                     organism=organism, pvalueCutoff=pvalue, limit = limit)
  kegg2 = list(enrichRes = kegg2,
               gridPlot = EnrichedView(EnrichedFilter(kegg2), top = 8, bottom = 0)
               + labs(title = "KEGG: Group2"))
  go2 = list(enrichRes = go2,
             gridPlot = EnrichedView(EnrichedFilter(go2), top = 8, bottom = 0)
             + labs(title = "Gene Ontology: Group2"))

  #====GO_KEGG_enrichment=====
  geneList = gg$diff[idx3]; names(geneList) = rownames(gg)[idx3]
  kegg3=EnrichAnalyzer(geneList = geneList, universe=universe,
                       method = enrich_method,type = "KEGG",
                       organism=organism, pvalueCutoff = pvalue, limit = limit)
  go3=EnrichAnalyzer(geneList = geneList, universe=universe,
                     method = enrich_method, type = "GOBP+GOMF+GOCC",
                     organism=organism, pvalueCutoff=pvalue, limit = limit)
  kegg3 = list(enrichRes = kegg3,
               gridPlot = EnrichedView(EnrichedFilter(kegg3), top = 0, bottom = 8)
               + labs(title = "KEGG: Group3"))
  go3 = list(enrichRes = go3,
             gridPlot = EnrichedView(EnrichedFilter(go3), top = 0, bottom = 8)
             + labs(title = "Gene Ontology: Group3"))

  #====GO_KEGG_enrichment=====
  geneList = gg$diff[idx4]; names(geneList) = rownames(gg)[idx4]
  kegg4=EnrichAnalyzer(geneList = geneList, universe=universe,
                       method = enrich_method, type = "KEGG",
                       organism=organism, pvalueCutoff = pvalue, limit = limit)
  go4=EnrichAnalyzer(geneList = geneList, universe=universe,
                     method = enrich_method, type = "GOBP+GOMF+GOCC",
                     organism=organism, pvalueCutoff=pvalue, limit = limit)
  kegg4 = list(enrichRes = kegg4,
               gridPlot = EnrichedView(EnrichedFilter(kegg4), top = 0, bottom = 8)
               + labs(title = "KEGG: Group4"))
  go4 = list(enrichRes = go4,
             gridPlot = EnrichedView(EnrichedFilter(go4), top = 0, bottom = 8)
             + labs(title = "Gene Ontology: Group4"))

  #====GO_KEGG_enrichment=====
  geneList = abs(gg$diff[idx1|idx2]); names(geneList) = rownames(gg)[idx1|idx2]
  kegg12=EnrichAnalyzer(geneList = geneList, universe=universe,
                        method = enrich_method, type = "KEGG",
                        organism=organism, pvalueCutoff = pvalue, limit = limit)
  go12=EnrichAnalyzer(geneList = geneList, universe=universe,
                      method = enrich_method, type = "GOBP+GOMF+GOCC",
                      organism=organism, pvalueCutoff = pvalue, limit = limit)
  kegg12 = list(enrichRes = kegg12,
                gridPlot = EnrichedView(EnrichedFilter(kegg12), top = 8, bottom = 0)
                + labs(title = "KEGG: Group1&2"))
  go12 = list(enrichRes = go12,
              gridPlot = EnrichedView(EnrichedFilter(go12), top = 8, bottom = 0)
              + labs(title = "Gene Ontology: Group1&2"))

  #====GO_KEGG_enrichment=====
  geneList = gg$diff[idx1|idx3]; names(geneList) = rownames(gg)[idx1|idx3]
  kegg13=EnrichAnalyzer(geneList = geneList, universe=universe,
                        method = enrich_method, type = "KEGG",
                        organism=organism, pvalueCutoff = pvalue, limit = limit)
  go13=EnrichAnalyzer(geneList = geneList, universe=universe,
                      method = enrich_method, type = "GOBP+GOMF+GOCC",
                      organism=organism, pvalueCutoff = pvalue, limit = limit)
  kegg13 = list(enrichRes = kegg13,
                gridPlot = EnrichedView(EnrichedFilter(kegg13), top = 4, bottom = 4)
               + labs(title = "KEGG: Group1&3"))
  go13 = list(enrichRes = go13,
              gridPlot = EnrichedView(EnrichedFilter(go13), top = 4, bottom = 4)
             + labs(title = "Gene Ontology: Group1&3"))

  #====GO_KEGG_enrichment=====
  geneList = gg$diff[idx4|idx2]; names(geneList) = rownames(gg)[idx4|idx2]
  kegg24=EnrichAnalyzer(geneList = geneList, universe=universe,
                        method = enrich_method,type = "KEGG",
                        organism=organism, pvalueCutoff = pvalue,
                        limit = limit)
  go24=EnrichAnalyzer(geneList = geneList, universe=universe,
                      method = enrich_method, type = "GOBP+GOMF+GOCC",
                      organism=organism, limit = limit, pvalueCutoff = pvalue)
  kegg24 = list(enrichRes = kegg24,
                gridPlot = EnrichedView(EnrichedFilter(kegg24), top = 4, bottom = 4)
                + labs(title = "KEGG: Group2&4"))
  go24 = list(enrichRes = go24,
              gridPlot = EnrichedView(EnrichedFilter(go24), top = 4, bottom = 4)
              + labs(title = "Gene Ontology: Group2&4"))

  #====GO_KEGG_enrichment=====
  geneList = gg$diff[idx3|idx4]; names(geneList) = rownames(gg)[idx3|idx4]
  kegg34=EnrichAnalyzer(geneList = geneList, universe=universe,
                        method = enrich_method, type = "KEGG", limit = limit,
                        organism=organism, pvalueCutoff = pvalue)
  go34=EnrichAnalyzer(geneList = geneList, universe=universe,
                      method = enrich_method, type = "GOBP+GOMF+GOCC",
                      organism=organism, limit = limit, pvalueCutoff = pvalue)
  kegg34 = list(enrichRes = kegg34,
                gridPlot = EnrichedView(EnrichedFilter(kegg34), top = 0, bottom = 8)
                + labs(title = "KEGG: Group3&4"))
  go34 = list(enrichRes = go34,
              gridPlot = EnrichedView(EnrichedFilter(go34), top = 0, bottom = 8)
              + labs(title = "Gene Ontology: Group3&4"))


  ###========Output results==================
  if(!is.null(filename)){
    if(!is.null(kegg1$enrichRes)){
      write.table(kegg1$enrichRes@result,
                  file.path(out.dir,paste0("Group1_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(kegg1$gridPlot,
             filename=file.path(out.dir,paste0("Group1_kegg_",filename,".png")),
             units = "in", width=width, height=height, ...)
    }
    if(!is.null(go1$enrichRes)){
      write.table(go1$enrichRes@result,
                  file.path(out.dir,paste0("Group1_go_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(go1$gridPlot,
             filename=file.path(out.dir,paste0("Group1_go_",filename,".png")),
             units = "in", width=width, height=height, ...)
    }
    if(!is.null(kegg2$enrichRes)){
      write.table(kegg2$enrichRes@result,
                  file.path(out.dir,paste0("Group2_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(kegg2$gridPlot,
             filename=file.path(out.dir,paste0("Group2_kegg_",filename,".png")),
             units = "in", width=width, height=height, ...)
    }
    if(!is.null(go2$enrichRes)){
        write.table(go2$enrichRes@result,
                    file.path(out.dir,paste0("Group2_go_",filename,".txt")),
                    sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
        ggsave(go2$gridPlot,
               filename=file.path(out.dir,paste0("Group2_go_",filename,".png")),
               units = "in", width=width, height=height, ...)
    }
    if(!is.null(kegg3$enrichRes)){
      write.table(kegg3$enrichRes@result,
                  file.path(out.dir,paste0("Group3_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(kegg3$gridPlot,
             filename=file.path(out.dir,paste0("Group3_kegg_",filename,".png")),
             units = "in", width=width, height=height, ...)
    }
    if(!is.null(go3$enrichRes)){
      write.table(go3$enrichRes@result,
                  file.path(out.dir,paste0("Group3_go_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(go3$gridPlot,
             filename=file.path(out.dir,paste0("Group3_go_",filename,".png")),
             units = "in", width=width, height=height, ...)
    }
    if(!is.null(kegg4$enrichRes)){
      write.table(kegg4$enrichRes@result,
                  file.path(out.dir,paste0("Group4_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(kegg4$gridPlot,
             filename=file.path(out.dir,paste0("Group4_kegg_",filename,".png")),
             units = "in", width=width, height=height, ...)
    }
    if(!is.null(go4$enrichRes)){
      write.table(go4$enrichRes@result,
                  file.path(out.dir,paste0("Group4_go_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(go4$gridPlot,
             filename=file.path(out.dir,paste0("Group4_go_",filename,".png")),
             units = "in", width=width, height=height, ...)
    }
    if(!is.null(kegg13$enrichRes)){
      write.table(kegg13$enrichRes@result,
                  file.path(out.dir,paste0("Group1&3_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(kegg13$gridPlot,
             filename=file.path(out.dir,paste0("Group1&3_kegg_",filename,".png")),
             units = "in", width=width, height=height, ...)
    }
    if(!is.null(go13$enrichRes)){
      write.table(go13$enrichRes@result,
                  file.path(out.dir,paste0("Group1&3_go_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(go13$gridPlot,
             filename=file.path(out.dir,paste0("Group1&3_go_",filename,".png")),
             units = "in", width=width, height=height, ...)
    }
    if(!is.null(kegg12$enrichRes)){
      write.table(kegg12$enrichRes@result,
                  file.path(out.dir,paste0("Group1&2_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(kegg12$gridPlot,
             filename=file.path(out.dir,paste0("Group1&2_kegg_",filename,".png")),
             units = "in", width=width, height=height, ...)
    }
    if(!is.null(go12$enrichRes)){
      write.table(go12$enrichRes@result,
                  file.path(out.dir,paste0("Group1&2_go_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(go12$gridPlot,
             filename=file.path(out.dir,paste0("Group1&2_go_",filename,".png")),
             units = "in", width=width, height=height, ...)
    }
    if(!is.null(kegg34$enrichRes)){
      write.table(kegg34$enrichRes@result,
                  file.path(out.dir,paste0("Group3&4_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(kegg34$gridPlot,
             filename=file.path(out.dir,paste0("Group3&4_kegg_",filename,".png")),
             units = "in", width=width, height=height, ...)
    }
    if(!is.null(go34$enrichRes)){
      write.table(go34$enrichRes@result,
                  file.path(out.dir,paste0("Group3&4_go_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(go34$gridPlot,
             filename=file.path(out.dir,paste0("Group3&4_go_",filename,".png")),
             units = "in", width=width, height=height, ...)
    }
    if(!is.null(kegg24$enrichRes)){
      write.table(kegg24$enrichRes@result,
                  file.path(out.dir,paste0("Group2&4_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(kegg24$gridPlot,
             filename=file.path(out.dir,paste0("Group2&4_kegg_",filename,".png")),
             units = "in", width=width, height=height, ...)
    }
    if(!is.null(go24$enrichRes)){
      write.table(go24$enrichRes@result,
                  file.path(out.dir,paste0("Group2&4_go_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(go24$gridPlot,
             filename=file.path(out.dir,paste0("Group2&4_go_",filename,".png")),
             units = "in", width=width, height=height, ...)
    }
  }
  return(list(kegg1=kegg1, kegg2=kegg2, kegg3=kegg3, kegg4=kegg4,
              kegg12=kegg12, kegg13=kegg13, kegg24=kegg24, kegg34=kegg34,
              go1=go1, go2=go2, go3=go3, go4=go4,
              go13=go13, go12=go12, go24=go24, go34=go34))
}

