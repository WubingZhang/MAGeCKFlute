#' Downstream analysis based on MAGeCK-RRA result
#'
#' Integrative analysis pipeline using the gene summary table in MAGeCK RRA results
#'
#' @docType methods
#' @name FluteRRA
#' @rdname FluteRRA
#' @aliases RRApipeline
#'
#' @param gene_summary a file path or a data frame, which has three columns named 'id', 'neg.fdr' and 'pos.fdr'.
#' @param prefix a character, indicating the prefix of output file name.
#' @param enrich_kegg One of "ORT"(Over-Representing Test), "GSEA"(Gene Set Enrichment Analysis), "DAVID",
#' "GOstats", and "HGT"(HyperGemetric test), or index from 1 to 5, specifying enrichment method used for kegg enrichment analysis.
#' @param organism a character, specifying organism, such as "hsa" or "Human"(default),
#' and "mmu" or "Mouse"
#' @param pvalueCutoff a numeric, specifying pvalue cutoff of enrichment analysis, default 1.
#' @param adjust one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param outdir output directory on disk
#'
#'
#' @author Wubing Zhang
#'
#' @note  See the vignette for an example of FluteRRA
#' The source can be found by typing \code{MAGeCKFlute:::FluteRRA}
#' or \code{getMethod("FluteRRA")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/FluteRRA.R}
#' Users should find it easy to customize this function.
#'
#' @return  All of the pipeline results is output into the \code{out.dir}/\code{prefix}_Results,
#' which includes a pdf file and a folder named 'RRA'.
#'
#' @details MAGeCK RRA allows for the comparison between two experimental conditions. It can identify
#' genes and sgRNAs are significantly selected between the two conditions. The most important output
#' of MAGeCK RRA is the file `gene_summary.txt`. MAGeCK RRA will output both the negative score and
#' positive score for each gene. A smaller score indicates higher gene importance. MAGeCK RRA  will
#' also output the statistical value for the scores of each gene. Genes that are significantly positively
#' and negatively selected can be identified based on the p-value or FDR.
#'
#' The downstream analysis of this function includes identifying positive and negative selection genes,
#' and performing biological functional category analysis and pathway enrichment analysis of these genes.
#'
#' @seealso \code{\link{FluteMLE}}
#'
#' @examples
#' \dontrun{
#' data(RRA_Data)
#' gene_summary = RRA_Data
#' # Run the FluteRRA pipeline
#' FluteRRA(gene_summary, prefix="BRAF", organism="hsa")
#' }
#'
#'
#' @export


#===read RRA results=====================================
FluteRRA <- function(gene_summary, prefix="Test", enrich_kegg="ORT",
                     organism="hsa", pvalueCutoff=0.25, adjust="BH",
                     outdir="."){
  #=========Prepare the running environment=========
  {
    loginfo("Create output dir and pdf file...")

    out.dir_sub=file.path(outdir, paste0(prefix, "_Flute_Results"))
    dir.create(file.path(out.dir_sub), showWarnings=FALSE)
    dir.create(file.path(out.dir_sub,"RRA"), showWarnings=FALSE)

    output_pdf = paste0(prefix,"_Flute.rra_summary.pdf")
    pdf(file.path(out.dir_sub, output_pdf),width=9,height = 4)
  }

  #=========Input data=========
  loginfo("Read RRA result ...")
  dd = ReadRRA(gene_summary, organism=organism)

  #enrichment analysis
  {
    universe=dd$ENTREZID
    idx=dd$neg.fdr<pvalueCutoff
    genes = dd[idx, "ENTREZID"]
    geneList=dd[idx, "neg.fdr"]
    names(geneList)=genes

    kegg.neg=enrichment_analysis(geneList=genes, universe=universe,
                                 method = enrich_kegg,type = "KEGG",
                                 organism=organism,pvalueCutoff=pvalueCutoff,
                                 plotTitle="KEGG: neg",color="#3f90f7",
                                 pAdjustMethod = adjust)
    bp.neg=enrichment_analysis(geneList=genes, universe=universe, method = "ORT",
                               type = "BP", organism=organism,
                               pvalueCutoff = pvalueCutoff, plotTitle="BP: neg",
                               color="#3f90f7", pAdjustMethod = adjust)

    grid.arrange(kegg.neg$gridPlot, bp.neg$gridPlot, ncol = 2)

    ggsave(kegg.neg$gridPlot,filename=
             file.path(out.dir_sub,"RRA/kegg.neg.png"),units = "in",
           width=400/100,height =270/100 )
    ggsave(bp.neg$gridPlot,filename=file.path(out.dir_sub,"RRA/bp.neg.png"),
           units = "in",width=400/100,height =270/100 )

    idx=dd$pos.fdr<pvalueCutoff
    genes = dd[idx, "ENTREZID"]
    geneList=dd[idx, "pos.fdr"]
    names(geneList)=genes

    kegg.pos=enrichment_analysis(geneList=genes, universe=universe,
                                 method = enrich_kegg, type = "KEGG",
                                 organism=organism, pvalueCutoff=pvalueCutoff,
                                 plotTitle="KEGG: pos",color="#e41a1c",
                                 pAdjustMethod = adjust)
    bp.pos=enrichment_analysis(geneList=genes, universe=universe, method = "ORT",
                               type = "BP", organism=organism,
                               pvalueCutoff = pvalueCutoff, plotTitle="BP: pos",
                               color="#e41a1c", pAdjustMethod = adjust)
    # gse=enrichment_analysis(geneList = geneList, genes=genes, method = "GSEA",
                              #type = "KEGG", pvalueCutoff = pvalueCutoff,
    #                         plotTitle="GSEA: RRA",color="#e41a1c",
    #                         pAdjustMethod = adjust)
    grid.arrange(kegg.pos$gridPlot, bp.pos$gridPlot, ncol = 2)

    ggsave(kegg.pos$gridPlot,filename=
             file.path(out.dir_sub,"RRA/kegg.pos.png"),units = "in",
           width=400/100,height =270/100 )
    ggsave(bp.pos$gridPlot,filename=
             file.path(out.dir_sub,"RRA/bp.pos.png"),units = "in",
           width=400/100,height =270/100 )

  }
  dev.off()
}
