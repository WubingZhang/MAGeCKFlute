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
#' \item{kegg13}{a list record enriched KEGG pathways for Group1&Group3 genes in 9-Square}
#' \item{kegg14}{a list record enriched KEGG pathways for Group1&Group4 genes in 9-Square}
#' \item{kegg23}{a list record enriched KEGG pathways for Group2&Group3 genes in 9-Square}
#' \item{kegg24}{a list record enriched KEGG pathways for Group2&Group4 genes in 9-Square}
#' \item{go1}{a list record enriched GO terms for Group1 genes in 9-Square}
#' \item{go2}{a list record enriched GO terms for Group2 genes in 9-Square}
#' \item{go3}{a list record enriched GO terms for Group3 genes in 9-Square}
#' \item{go4}{a list record enriched GO terms for Group4 genes in 9-Square}
#' \item{go13}{a list record enriched GO terms for Group1&Group3 genes in 9-Square}
#' \item{go14}{a list record enriched GO terms for Group1&Group4 genes in 9-Square}
#' \item{go23}{a list record enriched GO terms for Group2&Group3 genes in 9-Square}
#' \item{go24}{a list record enriched GO terms for Group2&Group4 genes in 9-Square}
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
#' dd = ReadBeta(mle.gene_summary, organism="hsa")
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

  idx=gg$group=="Group1"
  genes = rownames(gg)[idx]
  geneList = gg$diff[idx]
  names(geneList) = genes
  universe = rownames(gg)

  #====GO_KEGG_enrichment=====
  kegg1=enrichment_analysis(geneList = geneList, universe=universe,
                            method = enrich_method, type = "KEGG",
                            organism=organism, pvalueCutoff = pvalue,
                            main="KEGG: Group1", limit = limit)
  go1=enrichment_analysis(geneList = geneList, universe=universe,
                          method = "ORT", type = "GOBP+GOMF+GOCC", organism=organism,
                          pvalueCutoff = pvalue, main="Gene Ontology: Group1",
                          limit = limit)


  idx=gg$group=="Group2"
  genes=rownames(gg)[idx]
  geneList = gg$diff[idx]
  names(geneList) = genes
  #====GO_KEGG_enrichment=====
  kegg2=enrichment_analysis(geneList = geneList, universe=universe,
                            method = enrich_method, type = "KEGG",
                            organism=organism, pvalueCutoff=pvalue,
                            main="KEGG: Group2", limit = limit)
  go2=enrichment_analysis(geneList = geneList, universe=universe, method="ORT",
                          type = "GOBP+GOMF+GOCC",organism=organism, pvalueCutoff=pvalue,
                          main="Gene Ontology: Group2", limit = limit)

  idx=gg$group=="Group3"
  genes=rownames(gg)[idx]
  geneList = -gg$diff[idx]
  names(geneList) = genes
  #====GO_KEGG_enrichment=====
  kegg3=enrichment_analysis(geneList = geneList, universe=universe,
                            method = enrich_method,type = "KEGG",
                            organism=organism, pvalueCutoff = pvalue, limit = limit,
                            main="KEGG: Group3")
  go3=enrichment_analysis(geneList = geneList, universe=universe, method = "ORT",
                          type = "GOBP+GOMF+GOCC",organism=organism, pvalueCutoff=pvalue, limit = limit,
                          main="Gene Ontology: Group3")

  idx = gg$group=="Group4"
  genes = rownames(gg)[idx]
  geneList = -gg$diff[idx]
  names(geneList) = genes
  #====GO_KEGG_enrichment=====
  kegg4=enrichment_analysis(geneList = geneList, universe=universe,
                            method = enrich_method,type = "KEGG",
                            organism=organism, pvalueCutoff = pvalue, limit = limit,
                            main="KEGG: Group4")
  go4=enrichment_analysis(geneList = geneList, universe=universe, method="ORT",
                          type = "GOBP+GOMF+GOCC",organism=organism, pvalueCutoff=pvalue, limit = limit,
                          main="Gene Ontology: Group4")

  idx1=gg$group=="Group1"
  idx2=gg$group=="Group3"
  idx = idx1|idx2
  genes=rownames(gg)[idx]
  geneList = abs(gg$diff[idx])
  names(geneList) = genes
  #====GO_KEGG_enrichment=====
  kegg13=enrichment_analysis(geneList = geneList, universe=universe,
                             method = enrich_method,type = "KEGG",
                             organism=organism, pvalueCutoff = pvalue, limit = limit,
                             main="KEGG: Group1&3")
  go13=enrichment_analysis(geneList = geneList, universe=universe,
                           method = "ORT", type = "GOBP+GOMF+GOCC",organism=organism, limit = limit,
                           pvalueCutoff = pvalue, main="Gene Ontology: Group1&3")


  idx1=gg$group=="Group1"
  idx2=gg$group=="Group4"
  idx = idx1|idx2
  genes=rownames(gg)[idx]
  geneList = abs(gg$diff[idx])
  names(geneList) = genes
  #====GO_KEGG_enrichment=====
  kegg14=enrichment_analysis(geneList = geneList, universe=universe,
                             method = enrich_method,type = "KEGG", limit = limit,
                             organism=organism, pvalueCutoff = pvalue,
                             main="KEGG: Group1&4")
  go14=enrichment_analysis(geneList = geneList, universe=universe,
                           method = "ORT", type = "GOBP+GOMF+GOCC",organism=organism, limit = limit,
                           pvalueCutoff = pvalue, main="Gene Ontology: Group1&4")

  idx1=gg$group=="Group2"
  idx2=gg$group=="Group3"
  idx = idx1|idx2
  genes=rownames(gg)[idx]
  geneList = abs(gg$diff[idx])
  names(geneList) = genes
  #====GO_KEGG_enrichment=====
  kegg23=enrichment_analysis(geneList = geneList, universe=universe,
                             method = enrich_method,type = "KEGG", limit = limit,
                             organism=organism, pvalueCutoff = pvalue,
                             main="KEGG: Group2&3")
  go23=enrichment_analysis(geneList = geneList, universe=universe,
                           method = "ORT", type = "GOBP+GOMF+GOCC",organism=organism, limit = limit,
                           pvalueCutoff = pvalue, main="Gene Ontology: Group2&3")

  idx1=gg$group=="Group2"
  idx2=gg$group=="Group4"
  idx = idx1|idx2
  genes=rownames(gg)[idx]
  geneList = abs(gg$diff[idx])
  names(geneList) = genes
  #====GO_KEGG_enrichment=====
  kegg24=enrichment_analysis(geneList = geneList, universe=universe,
                             method = enrich_method,type = "KEGG",
                             organism=organism,pvalueCutoff = pvalue, limit = limit,
                             main="KEGG: Group2&4")
  go24=enrichment_analysis(geneList = geneList, universe=universe,
                           method = "ORT", type = "GOBP+GOMF+GOCC",organism=organism, limit = limit,
                           pvalueCutoff = pvalue, main="Gene Ontology: Group2&4")

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
    if(!is.null(kegg14$enrichRes)){
      write.table(kegg14$enrichRes@result,
                  file.path(out.dir,paste0("Group1&4_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(kegg14$gridPlot,
             filename=file.path(out.dir,paste0("Group1&4_kegg_",filename,".png")),
             units = "in", width=width, height=height, ...)
    }
    if(!is.null(go14$enrichRes)){
      write.table(go14$enrichRes@result,
                  file.path(out.dir,paste0("Group1&4_go_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(go14$gridPlot,
             filename=file.path(out.dir,paste0("Group1&4_go_",filename,".png")),
             units = "in", width=width, height=height, ...)
    }
    if(!is.null(kegg23$enrichRes)){
      write.table(kegg23$enrichRes@result,
                  file.path(out.dir,paste0("Group2&3_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(kegg23$gridPlot,
             filename=file.path(out.dir,paste0("Group2&3_kegg_",filename,".png")),
             units = "in", width=width, height=height, ...)
    }
    if(!is.null(go23$enrichRes)){
      write.table(go23$enrichRes@result,
                  file.path(out.dir,paste0("Group2&3_go_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(go23$gridPlot,
             filename=file.path(out.dir,paste0("Group2&3_go_",filename,".png")),
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
              kegg13=kegg13, kegg14=kegg14, kegg23=kegg23, kegg24=kegg24,
              go1=go1, go2=go2, go3=go3, go4=go4,
              go13=go13, go14=go14, go23=go23, go24=go24))
}

