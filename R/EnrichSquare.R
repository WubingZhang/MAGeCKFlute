#' Enrichment analysis for selected treatment related genes
#'
#' Do enrichment analysis for selected treatment related genes in 9-squares
#'
#' @docType methods
#' @name EnrichSquare
#' @rdname EnrichSquare
#'
#' @param beta Data frame, which has columns of 'ENTREZID' and 'group'.
#' @param pvalue Pvalue cutoff.
#' @param enrich_method One of "ORT"(Over-Representing Test), "DAVID", "GOstats", and "HGT"(HyperGemetric test).
#' @param organism A character, specifying organism, such as "hsa" or "Human"(default), and "mmu" or "Mouse"
#' @param adjust One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param filename Suffix of output file name. NULL(default) means no output.
#' @param out.dir Path to save plot to (combined with filename).
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
#' \item{kegg1234}{a list record enriched KEGG pathways for Group1&Group2&Group3&Group4 genes in 9-Square}
#' \item{bp1}{a list record enriched GO BP terms for Group1 genes in 9-Square}
#' \item{bp2}{a list record enriched GO BP terms for Group2 genes in 9-Square}
#' \item{bp3}{a list record enriched GO BP terms for Group3 genes in 9-Square}
#' \item{bp4}{a list record enriched GO BP terms for Group4 genes in 9-Square}
#' \item{bp13}{a list record enriched GO BP terms for Group1&Group3 genes in 9-Square}
#' \item{bp14}{a list record enriched GO BP terms for Group1&Group4 genes in 9-Square}
#' \item{bp23}{a list record enriched GO BP terms for Group2&Group3 genes in 9-Square}
#' \item{bp24}{a list record enriched GO BP terms for Group2&Group4 genes in 9-Square}
#' \item{bp1234}{a list record enriched GO BP terms for Group1&Group2&Group3&Group4 genes in 9-Square}
#' Each item in the returned list has two sub items:
#' \item{gridPlot}{an object created by \code{ggplot}, which can be assigned and further customized.}
#' \item{enrichRes}{a enrichResult instance.}
#'
#' @author Wubing Zhang
#'
#' @note  See the vignette for an example of EnrichSquare
#' The source can be found by typing \code{MAGeCKFlute:::EnrichSquare}
#' or \code{getMethod("EnrichSquare")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/EnrichSquare.R}
#' Users should find it easy to customize this function.
#'
#' @seealso \code{\link{SquareView}}
#' @seealso \code{\link{EnrichSquare}}
#'
# @examples
# \dontrun{
#   data(MLE_Data)
#'  # Read beta score from gene summary table in MAGeCK MLE results
#   dd = ReadBeta(MLE_Data, organism="hsa")
#   E1 = EnrichSquare(dd, ctrlname = "D7_R1", treatname = "PLX7_R1", pvalue=0.05, adjust="BH", enrich_method="ORT", organism="hsa")
#   print(E1$kegg1$gridPlot)
# }
#


#enrichment for square grouped genes
EnrichSquare <- function(beta, pvalue = 1,enrich_method="ORT",
                         organism="hsa", adjust="BH", filename=NULL,
                         out.dir="."){
  loginfo("Enrichment analysis of 9 Square grouped genes ...")

  ## ===========Enrichment===================
  gg=beta

  idx=gg$group=="Group1"
  genes=as.character(gg$ENTREZID[idx])
  universe=as.character(gg$ENTREZID)

  #====GO_KEGG_enrichment=====
  kegg1=enrichment_analysis(geneList = genes, universe=universe,
                            method = enrich_method,type = "KEGG",
                            organism=organism, pvalueCutoff = pvalue,color="#e41a1c",
                            plotTitle="KEGG: Group1", pAdjustMethod = adjust)
  bp1=enrichment_analysis(geneList = genes, universe=universe,
                          method = "ORT", type = "BP",organism=organism,
                          pvalueCutoff = pvalue, plotTitle="BP: Group1",
                          color="#e41a1c", pAdjustMethod = adjust)


  idx=gg$group=="Group2"
  genes=as.character(gg$ENTREZID[idx])
  #====GO_KEGG_enrichment=====
  kegg2=enrichment_analysis(geneList = genes, universe=universe,
                            method = enrich_method,type = "KEGG",
                            organism=organism, pvalueCutoff=pvalue, color="#e41a1c",
                            plotTitle="KEGG: Group2", pAdjustMethod=adjust)
  bp2=enrichment_analysis(geneList = genes, universe=universe, method="ORT",
                          type = "BP",organism=organism, pvalueCutoff=pvalue, color="#e41a1c",
                          plotTitle="BP: Group2", pAdjustMethod = adjust)

  idx=gg$group=="Group3"
  genes=as.character(gg$ENTREZID[idx])
  #====GO_KEGG_enrichment=====
  kegg3=enrichment_analysis(geneList = genes, universe=universe,
                            method = enrich_method,type = "KEGG",
                            organism=organism, pvalueCutoff = pvalue,
                            color="#3f90f7", plotTitle="KEGG: Group3", pAdjustMethod = adjust)
  bp3=enrichment_analysis(geneList = genes, universe=universe, method = "ORT",
                          type = "BP",organism=organism, pvalueCutoff=pvalue,
                          color="#3f90f7", plotTitle="BP: Group3", pAdjustMethod = adjust)

  idx=gg$group=="Group4"
  genes=as.character(gg$ENTREZID[idx])
  #====GO_KEGG_enrichment=====
  kegg4=enrichment_analysis(geneList = genes, universe=universe,
                            method = enrich_method,type = "KEGG",
                            organism=organism, pvalueCutoff = pvalue,
                            color="#3f90f7", plotTitle="KEGG: Group4", pAdjustMethod=adjust)
  bp4=enrichment_analysis(geneList = genes, universe=universe, method="ORT",
                          type = "BP",organism=organism, pvalueCutoff=pvalue,
                          color="#3f90f7", plotTitle="BP: Group4", pAdjustMethod = adjust)

  idx1=gg$group=="Group1"
  idx2=gg$group=="Group3"
  idx = idx1|idx2
  genes=as.character(gg$ENTREZID[idx])
  #====GO_KEGG_enrichment=====
  kegg13=enrichment_analysis(geneList = genes, universe=universe,
                             method = enrich_method,type = "KEGG",
                             organism=organism, pvalueCutoff = pvalue,
                             color="#6daf61", plotTitle="KEGG: Group1&3", pAdjustMethod = adjust)
  bp13=enrichment_analysis(geneList = genes, universe=universe,
                           method = "ORT", type = "BP",organism=organism,
                           pvalueCutoff = pvalue, plotTitle="BP: Group1&3",
                           color="#6daf61", pAdjustMethod = adjust)


  idx1=gg$group=="Group1"
  idx2=gg$group=="Group4"
  idx = idx1|idx2
  genes=as.character(gg$ENTREZID[idx])
  #====GO_KEGG_enrichment=====
  kegg14=enrichment_analysis(geneList = genes, universe=universe,
                             method = enrich_method,type = "KEGG",
                             organism=organism, pvalueCutoff = pvalue,
                             color="#6daf61", plotTitle="KEGG: Group1&4", pAdjustMethod=adjust)
  bp14=enrichment_analysis(geneList = genes, universe=universe,
                           method = "ORT", type = "BP",organism=organism,
                           pvalueCutoff = pvalue, plotTitle="BP: Group1&4",
                           color="#6daf61", pAdjustMethod = adjust)

  idx1=gg$group=="Group2"
  idx2=gg$group=="Group3"
  idx = idx1|idx2
  genes=as.character(gg$ENTREZID[idx])
  #====GO_KEGG_enrichment=====
  kegg23=enrichment_analysis(geneList = genes, universe=universe,
                             method = enrich_method,type = "KEGG",
                             organism=organism, pvalueCutoff = pvalue,
                             color="#6daf61", plotTitle="KEGG: Group2&3", pAdjustMethod = adjust)
  bp23=enrichment_analysis(geneList = genes, universe=universe,
                           method = "ORT", type = "BP",organism=organism,
                           pvalueCutoff = pvalue, plotTitle="BP: Group2&3",
                           color="#6daf61", pAdjustMethod = adjust)

  idx1=gg$group=="Group2"
  idx2=gg$group=="Group4"
  idx = idx1|idx2
  genes=as.character(gg$ENTREZID[idx])

  #====GO_KEGG_enrichment=====
  kegg24=enrichment_analysis(geneList = genes, universe=universe,
                             method = enrich_method,type = "KEGG",
                             organism=organism,pvalueCutoff = pvalue,
                             color="#6daf61", plotTitle="KEGG: Group2&4", pAdjustMethod = adjust)
  bp24=enrichment_analysis(geneList = genes, universe=universe,
                           method = "ORT", type = "BP",organism=organism,
                           pvalueCutoff = pvalue, plotTitle="BP: Group2&4",
                           color="#6daf61", pAdjustMethod = adjust)

  idx1=gg$group=="Group1"
  idx2=gg$group=="Group2"
  idx3=gg$group=="Group3"
  idx4=gg$group=="Group4"
  idx = idx1|idx2|idx3|idx4
  genes=as.character(gg$ENTREZID[idx])
  #====GO_KEGG_enrichment=====
  kegg1234=enrichment_analysis(geneList = genes, universe=universe,
                               method = enrich_method,type = "KEGG",
                               organism=organism, pvalueCutoff = pvalue,
                               plotTitle="KEGG: Group1&2&3&4",
                               color="#6daf61", pAdjustMethod = adjust)
  bp1234=enrichment_analysis(geneList = genes, universe=universe,
                             method = "ORT", type = "BP",organism=organism,
                             pvalueCutoff = pvalue, plotTitle="BP: Group1&2&3&4",
                             color="#6daf61", pAdjustMethod = adjust)
  ###========Output results=================================================
  if(!is.null(filename)){
    if(!is.null(kegg1$enrichRes)){
      write.table(kegg1$enrichRes@result,
                  file.path(out.dir,paste0("Group1_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(kegg1$gridPlot,
             filename=file.path(out.dir,paste0("Group1_kegg_",filename,".png")),
             units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(bp1$enrichRes)){
      write.table(bp1$enrichRes@result,
                  file.path(out.dir,paste0("Group1_bp_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(bp1$gridPlot,
             filename=file.path(out.dir,paste0("Group1_bp_",filename,".png")),
             units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(kegg2$enrichRes)){
      write.table(kegg2$enrichRes@result,
                  file.path(out.dir,paste0("Group2_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(kegg2$gridPlot,
             filename=file.path(out.dir,paste0("Group2_kegg_",filename,".png")),
             units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(bp2$enrichRes)){
        write.table(bp2$enrichRes@result,
                    file.path(out.dir,paste0("Group2_bp_",filename,".txt")),
                    sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
        ggsave(bp2$gridPlot,
               filename=file.path(out.dir,paste0("Group2_bp_",filename,".png")),
               units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(kegg3$enrichRes)){
      write.table(kegg3$enrichRes@result,
                  file.path(out.dir,paste0("Group3_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(kegg3$gridPlot,
             filename=file.path(out.dir,paste0("Group3_kegg_",filename,".png")),
             units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(bp3$enrichRes)){
      write.table(bp3$enrichRes@result,
                  file.path(out.dir,paste0("Group3_bp_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(bp3$gridPlot,
             filename=file.path(out.dir,paste0("Group3_bp_",filename,".png")),
             units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(kegg4$enrichRes)){
      write.table(kegg4$enrichRes@result,
                  file.path(out.dir,paste0("Group4_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(kegg4$gridPlot,
             filename=file.path(out.dir,paste0("Group4_kegg_",filename,".png")),
             units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(bp4$enrichRes)){
      write.table(bp4$enrichRes@result,
                  file.path(out.dir,paste0("Group4_bp_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(bp4$gridPlot,
             filename=file.path(out.dir,paste0("Group4_bp_",filename,".png")),
             units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(kegg13$enrichRes)){
      write.table(kegg13$enrichRes@result,
                  file.path(out.dir,paste0("Group1&3_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(kegg13$gridPlot,
             filename=file.path(out.dir,paste0("Group1&3_kegg_",filename,".png")),
             units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(bp13$enrichRes)){
      write.table(bp13$enrichRes@result,
                  file.path(out.dir,paste0("Group1&3_bp_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(bp13$gridPlot,
             filename=file.path(out.dir,paste0("Group1&3_bp_",filename,".png")),
             units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(kegg14$enrichRes)){
      write.table(kegg14$enrichRes@result,
                  file.path(out.dir,paste0("Group1&4_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(kegg14$gridPlot,
             filename=file.path(out.dir,paste0("Group1&4_kegg_",filename,".png")),
             units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(bp14$enrichRes)){
      write.table(bp14$enrichRes@result,
                  file.path(out.dir,paste0("Group1&4_bp_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(bp14$gridPlot,
             filename=file.path(out.dir,paste0("Group1&4_bp_",filename,".png")),
             units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(kegg23$enrichRes)){
      write.table(kegg23$enrichRes@result,
                  file.path(out.dir,paste0("Group2&3_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(kegg23$gridPlot,
             filename=file.path(out.dir,paste0("Group2&3_kegg_",filename,".png")),
             units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(bp23$enrichRes)){
      write.table(bp23$enrichRes@result,
                  file.path(out.dir,paste0("Group2&3_bp_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(bp23$gridPlot,
             filename=file.path(out.dir,paste0("Group2&3_bp_",filename,".png")),
             units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(kegg24$enrichRes)){
      write.table(kegg24$enrichRes@result,
                  file.path(out.dir,paste0("Group2&4_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(kegg24$gridPlot,
             filename=file.path(out.dir,paste0("Group2&4_kegg_",filename,".png")),
             units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(bp24$enrichRes)){
      write.table(bp24$enrichRes@result,
                  file.path(out.dir,paste0("Group2&4_bp_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(bp24$gridPlot,
             filename=file.path(out.dir,paste0("Group2&4_bp_",filename,".png")),
             units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(kegg1234$enrichRes)){
      write.table(kegg1234$enrichRes@result,
                  file.path(out.dir,paste0("Group1&2&3&4_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(kegg1234$gridPlot,
             filename=file.path(out.dir,paste0("Group1&2&3&4_kegg_",filename,".png")),
             units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(bp1234$enrichRes)){
      write.table(bp1234$enrichRes@result,
                  file.path(out.dir,paste0("Group1&2&3&4_bp_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(bp1234$gridPlot,
             filename=file.path(out.dir,paste0("Group1&2&3&4_bp_",filename,".png")),
             units = "in",width=400/100,height =270/100 )
    }
  }
  return(list(kegg1=kegg1, kegg2=kegg2, kegg3=kegg3, kegg4=kegg4,
              kegg13=kegg13, kegg14=kegg14, kegg23=kegg23, kegg24=kegg24,
              kegg1234=kegg1234, bp1=bp1, bp2=bp2, bp3=bp3, bp4=bp4,
              bp13=bp13, bp14=bp14, bp23=bp23, bp24=bp24, bp1234=bp1234))

}

