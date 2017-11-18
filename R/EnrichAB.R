#' Enrichment analysis for Positive and Negative selection genes
#'
#' Do enrichment analysis for selected genes, in which positive selection and negative selection
#' are termed as GroupA and GroupB
#'
#' @docType methods
#' @name enrichAB
#' @rdname enrichAB
#'
#' @param data a data frame containing columns of "ENTREZID" and "diff".
#' @param pvalue pvalue cutoff.
#' @param enrich_method One of "ORT"(Over-Representing Test), "GSEA"(Gene Set Enrichment Analysis), "DAVID",
#' "GOstats", and "HGT"(HyperGemetric test), or index from 1 to 5
#' @param organism a character, specifying organism, such as "hsa" or "Human"(default), and "mmu" or "Mouse"

#' @param adjust one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param filename suffix of output file name. NULL(default) means no output.
#' @param out.dir Path to save plot to (combined with filename).
#' @param gsea boolean, specifying if do GSEA for GroupA and GroupB genes. Default gsea=FALSE.
#'
#' @return a list containing enrichment results for each group genes. This list contains items four
#' items, \code{keggA}, \code{keggB}, \code{bpA}, \code{bpB}. Four items are all list object, containing
#' subitems of \code{gridPlot} and \code{enrichRes}. \code{gridPlot} is a ggplot object, and
#' \code{enrichRes} is a enrichResult instance
#'
#' @author Binbin Wang
#'
#' @note  See the vignette for an example of EnrichAB
#' The source can be found by typing \code{MAGeCKFlute:::EnrichAB}
#' or \code{getMethod("EnrichAB")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/EnrichAB.R}
#' Users should find it easy to customize this function.
#'
#' @seealso \code{\link{EnrichSquare}}
#'
# @examples
#' data(MLE_Data)
#' # Read beta score from gene summary table in MAGeCK MLE results
# dd = ReadBeta(MLE_Data, ctrlName = "D7_R1", treatName = "PLX7_R1", organism="hsa")
#
# \dontrun{
#   data=ScatterView(dd)$data
#   #BP and KEGG enrichment analysis
#   enrich_result = EnrichAB(data, pvalue=0.05, organism="hsa")
#   print(enrich_result$keggA$gridPlot)
#   print(enrich_result$bpA$gridPlot)
# }
#

# enrichment for GroupA and GrouB genes
EnrichAB <- function(data, pvalue=0.05, enrich_method="ORT",
                     organism="hsa", adjust="BH", filename=NULL,
                     out.dir=".", gsea=FALSE){
  loginfo("Enrichment analysis of GroupA and GroupB genes ...")
  gg=data
  ##===================enrichment for GroupA==============================
  idx1=gg$group=="up"
  genes=as.character(gg$ENTREZID[idx1])
  geneList=gg$diff[idx1]
  names(geneList)=genes
  universe=as.character(gg$ENTREZID)

  #====GO_KEGG_enrichment=====
  keggA=enrichment_analysis(geneList = genes, universe=universe,
                            method = enrich_method,type = "KEGG",
                            organism=organism,pvalueCutoff = pvalue,
                            plotTitle="KEGG: GroupA",gridColour="#e41a1c",
                            pAdjustMethod = adjust)
  bpA=enrichment_analysis(geneList = genes, universe=universe,
                          method = "ORT", type = "BP", organism=organism,
                          pvalueCutoff = pvalue, plotTitle="BP: GroupA",
                          gridColour="#e41a1c", pAdjustMethod = adjust)
  if(gsea){
    requireNamespace("clusterProfiler", quietly=TRUE) || stop("need clusterProfiler package")
    gseA=enrichment_analysis(geneList = geneList, method = "GSEA",
                             type = "KEGG", organism=organism,
                             pvalueCutoff = pvalue, plotTitle="GSEA: GroupA",
                             gridColour="#e41a1c", pAdjustMethod = adjust)
  }
  ##=============Enrichment for GroupB========================================
  idx2=gg$group=="down"
  genes=gg$ENTREZID[idx2]
  geneList=gg$diff
  names(geneList)=gg$ENTREZID
  #====GO_KEGG_enrichment=====
  keggB=enrichment_analysis(geneList = genes, universe=universe,
                            method = enrich_method,type = "KEGG",
                            organism=organism, pvalueCutoff = pvalue,
                            plotTitle="KEGG: GroupB",gridColour="#377eb8",
                            pAdjustMethod = adjust)
  bpB = enrichment_analysis(geneList = genes, universe=universe,
                            method = "ORT",type = "BP",organism=organism,
                            pvalueCutoff = pvalue, plotTitle="BP: GroupB",
                            gridColour="#377eb8", pAdjustMethod = adjust)
  if(gsea){
    gseB=enrichment_analysis(geneList = geneList, method = "GSEA",
                             type = "KEGG", organism=organism,
                             pvalueCutoff = pvalue, plotTitle="GSEA: GroupB",
                             gridColour="#377eb8", pAdjustMethod = adjust)
  }
  ##================output results=============================================
  if(!is.null(filename)){
    ####===========GSEA results===================================
    if(gsea){
      p1=ggplot()
      p1=p1+geom_text(aes(x=0,y=0,label="No enriched terms"),size=6)
      p1=p1+theme_void()

      dir.create(file.path(out.dir,"GSEA_results"), showWarnings=FALSE)
      ##=========GroupA GSEA plot=================================
      if(!is.null(gseA$enrichRes) && nrow(gseA$enrichRes@result)>0){
        for(term in gseA$enrichRes@result$ID[nrow(gseA$enrichRes@result):1]){
          png(file.path(out.dir,paste0("GSEA_results/GroupA_gse_",
                                       term, "_", filename,".png")),
              units = "in",width=400/100,height =270/100,res=300)
          p1 = gseaplot(gseA$enrichRes, term)$runningScore
          dev.off()
        }
        p1 <- p1+xlab("Ranked list of genes")+ylab("Enrichment score")
        p1 <- p1+labs(title=as.character(gseA$enrichRes@result$Description[1]))
        p1 <- p1+theme(axis.text.x=element_text(size=6, face="plain",
                                                colour='black'))
        p1 <- p1+theme(axis.text.y=element_text(size=6, face="plain",
                                                colour='black'))
        p1=p1+theme(plot.title = element_text(hjust = 0.5,size=10,
                                              face="plain", colour='black'))
        p1 <- p1+theme(panel.grid.minor=element_blank(),
                       panel.background=element_blank())
        write.table(gseA$enrichRes@result,
              file.path(out.dir,
                  paste0("GSEA_results/GroupA_gse_",filename,".txt")),
              sep="\t", row.names = FALSE,col.names = TRUE,quote= FALSE)
      }
      gseA$gseaplot = p1
      ##=========GroupB GSEA plot==================================
      if(!is.null(gseB$enrichRes) && nrow(gseB$enrichRes@result)>0){
        for(term in gseB$enrichRes@result$ID[nrow(gseB$enrichRes@result):1]){
          png(file.path(out.dir,paste0("GSEA_results/GroupB_gse_",
                                       term, "_", filename,".png")),
              units = "in",width=400/100,height =270/100,res=300)
          p1 = gseaplot(gseB$enrichRes, term)$runningScore
          dev.off()
        }
        p1 <- p1+xlab("Ranked list of genes")+ylab("Enrichment score")
        p1 <- p1+labs(title=as.character(gseB$enrichRes@result$Description[1]))
        p1 <- p1+theme(axis.text.x=element_text(size=6, face="plain",
                                                colour='black'))
        p1 <- p1+theme(axis.text.y=element_text(size=6, face="plain",
                                                colour='black'))
        p1=p1+theme(plot.title = element_text(hjust = 0.5,size=10,
                                              face="plain", colour='black'))
        p1 <- p1+theme(panel.grid.minor=element_blank(),
                       panel.background=element_blank())
        write.table(gseB$enrichRes@result,
                    file.path(out.dir,paste0("GSEA_results/GroupB_gse_",
                                             filename,".txt")),
                    sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      }
      gseB$gseaplot = p1
    }
    ##=========Save GroupA enrichment results===========================
    if(!is.null(keggA$enrichRes)){
      write.table(keggA$enrichRes@result,
                  file.path(out.dir,paste0("GroupA_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(keggA$gridPlot,
             filename=file.path(out.dir,paste0("GroupA_kegg_",
                                               filename,".png")),
             units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(bpA$enrichRes)){
      write.table(bpA$enrichRes@result,
                  file.path(out.dir,paste0("GroupA_bp_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(bpA$gridPlot,
             filename=file.path(out.dir,paste0("GroupA_bp_",filename,".png")),
             units = "in",width=400/100,height =270/100 )
    }
    ##=========Save GroupB enrichment results===========================
    if(!is.null(keggB$enrichRes)){
      write.table(keggB$enrichRes@result,
                  file.path(out.dir,paste0("GroupB_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(keggB$gridPlot,
             filename=file.path(out.dir,paste0("GroupB_kegg_",filename,".png")),
             units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(bpB$enrichRes)){
      write.table(bpB$enrichRes@result,
                  file.path(out.dir,paste0("GroupB_bp_",filename,".txt")),
                  sep="\t", row.names = FALSE,col.names = TRUE,quote=FALSE)
      ggsave(bpB$gridPlot,
             filename=file.path(out.dir,paste0("GroupB_bp_",filename,".png")),
             units = "in",width=400/100,height =270/100 )
    }
  }
  ##=========Return results=====================================
  if(gsea){
    return(list(keggA=keggA, bpA=bpA, gseA=gseA,
                keggB=keggB, bpB=bpB, gseB=gseB))
  }else{
    return(list(keggA=keggA, bpA=bpA, keggB=keggB, bpB=bpB))
  }
}

