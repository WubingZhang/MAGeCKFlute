#' Downstream analysis based on MAGeCK-RRA result
#'
#' Integrative analysis pipeline using the gene summary table in MAGeCK RRA results
#'
#' @docType methods
#' @name FluteRRA
#' @rdname FluteRRA
#' @aliases RRApipeline
#'
#' @param gene_summary A file path or a data frame of gene summary data.
#' @param sgrna_summary A file path or a data frame of sgRNA summary data.
#' @param lfcCutoff A two-length vector (default: c(-1, 1)), specifying the logFC cutoff
#' for negative selection and positive selection.
#' @param organism "hsa" or "mmu".
#' @param pathway_limit A two-length vector (default: c(3, 50)), specifying the minimal and
#' maximal size of gene sets for enrichent analysis.
#' @param pvalueCutoff A numeric, specifying pvalue cutoff of enrichment analysis, default 1.
#' @param adjust One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param prefix A character, indicating the prefix of output file name.
#' @param width The width of summary pdf in inches.
#' @param height The height of summary pdf in inches.
#' @param outdir Output directory on disk.
#'
#' @author Wubing Zhang
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
#' data(rra.gene_summary)
#' \dontrun{
#'     # Run the FluteRRA pipeline
#'     FluteRRA(rra.gene_summary, prefix="GSC", organism="hsa")
#' }
#'
#'
#' @export

#===read RRA results=====================================
FluteRRA <- function(gene_summary, sgrna_summary, lfcCutoff = c(-1, 1), organism = "hsa",
                     pathway_limit = c(3, 50), pvalueCutoff = 0.25, adjust = "BH",
                     prefix = "Test", width = 12, height = 6, outdir = "."){
  #=========Prepare the running environment=========
  {
    message(Sys.time(), " # Create output dir and pdf file ...")
    prefix = gsub(".*\\/|.*\\\\", "", prefix)
    out.dir_sub=file.path(outdir, paste0(prefix, "_Flute_Results"))
    dir.create(file.path(out.dir_sub), showWarnings=FALSE)
    dir.create(file.path(out.dir_sub,"RRA"), showWarnings=FALSE)

    output_pdf = paste0(prefix,"_Flute.rra_summary.pdf")
    pdf(file.path(out.dir_sub, output_pdf),width=width, height = height)
  }

  #=========Input data=========
  message(Sys.time(), " # Read RRA result ...")
  dd = ReadRRA(gene_summary, organism=organism)
  dd.sgrna = ReadsgRRA(sgrna_summary)
  # Volcano plot
  p1 = VolcanoView(dd, x = "LFC", y = "FDR", Label = "Official")
  ggsave(p1, filename = file.path(out.dir_sub,"RRA/VolcanoView_RRA.png"),
         units = "in", width = 6.5, height = 4)
  geneList= dd$LFC
  names(geneList) = dd$Official
  p3 = RankView(geneList)
  ggsave(p3, filename = file.path(out.dir_sub,"RRA/RankView_Gene.png"),
         units = "in", width = 6.5, height = 4)
  topgenes = as.vector(p1$data$Label[!is.na(p1$data$Label) & p1$data$LFC>0])
  botgenes = as.vector(p1$data$Label[!is.na(p1$data$Label) & p1$data$LFC<0])
  p2 = sgRankView(dd.sgrna, top = 0, bottom = 0, gene = topgenes)
  ggsave(p2, filename = file.path(out.dir_sub,"RRA/RankView_sgRNA_top8_.png"),
         units = "in", width = 6.5, height = 4)
  p4 = sgRankView(dd.sgrna, top = 0, bottom = 0, gene = botgenes)
  ggsave(p4, filename = file.path(out.dir_sub,"RRA/RankView_sgRNA_bott8_.png"),
         units = "in", width = 6.5, height = 4)
  grid.arrange(p1, p3, p2, p4, ncol = 2)
  #enrichment analysis
  {
    universe=dd$EntrezID
    geneList= dd$LFC
    names(geneList) = dd$EntrezID
    idx=dd$LFC<lfcCutoff[1]
    kegg.neg=enrichment_analysis(geneList=geneList[idx], universe=universe,
                                 method = "HGT", type = "KEGG+BIOCARTA+REACTOME+EHMN+PID+WikiPathways",
                                 organism=organism,pvalueCutoff=pvalueCutoff,
                                 plotTitle="Pathway: neg",color="#3f90f7",
                                 pAdjustMethod = adjust, limit = pathway_limit)
    go.neg=enrichment_analysis(geneList=geneList[idx], universe=universe, method = "HGT",
                               type = "GOBP+GPMF+GOCC", organism=organism,
                               pvalueCutoff = pvalueCutoff, plotTitle="Gene Ontology: neg",
                               color="#3f90f7", pAdjustMethod = adjust, limit = pathway_limit)

    ggsave(filename = file.path(out.dir_sub, "RRA/kegg.neg.png"), kegg.neg$gridPlot,
           units = "in", width = 6.5, height = 4)
    saveRDS(kegg.neg, file.path(out.dir_sub, "RRA/EnrichRes_kegg.neg.rds"))
    ggsave(filename=file.path(out.dir_sub, "RRA/go.neg.png"), go.neg$gridPlot,
           units = "in", width = 6.5, height = 4)
    saveRDS(go.neg, file.path(out.dir_sub, "RRA/EnrichRes_go.neg.rds"))

    idx=dd$LFC>lfcCutoff[2]
    kegg.pos=enrichment_analysis(geneList=geneList[idx], universe=universe,
                                 method = "HGT", type = "KEGG+BIOCARTA+REACTOME+EHMN+PID+WikiPathways",
                                 organism=organism, pvalueCutoff=pvalueCutoff,
                                 plotTitle="Pathway: pos",color="#e41a1c",
                                 pAdjustMethod = adjust, limit = pathway_limit)
    go.pos=enrichment_analysis(geneList=geneList[idx], universe=universe, method = "HGT",
                               type = "GOBP+GOMF+GOCC", organism=organism,
                               pvalueCutoff = pvalueCutoff, plotTitle="Gene Ontology: pos",
                               color="#e41a1c", pAdjustMethod = adjust, limit = pathway_limit)
    ggsave(kegg.pos$gridPlot,filename=file.path(out.dir_sub,"RRA/kegg.pos.png"),
           units = "in", width = 6.5, height = 4)
    saveRDS(kegg.pos, file.path(out.dir_sub, "RRA/EnrichRes_kegg.pos.rds"))
    ggsave(go.pos$gridPlot,filename=file.path(out.dir_sub,"RRA/go.pos.png"),
           units = "in", width = 6.5, height = 4)
    saveRDS(go.pos, file.path(out.dir_sub, "RRA/EnrichRes_go.pos.rds"))

    grid.arrange(kegg.neg$gridPlot, kegg.pos$gridPlot, go.neg$gridPlot, go.pos$gridPlot, ncol = 2)

    if(!is.null(kegg.neg$enrichRes)){
      p1 = EnrichedGeneView(enrichment=kegg.neg$enrichRes@result, geneList, keytype = "Entrez",
                            gene_cutoff = lfcCutoff, top = 0, bottom = 15)
      ggsave(filename = file.path(out.dir_sub,"RRA/EnrichedGeneView_kegg.neg.png"), p1,
             units = "in", width = 10, height = 7)
      grid.arrange(p1, ncol = 1)
    }
    if(!is.null(go.neg$enrichRes)){
      p2 = EnrichedGeneView(enrichment=go.neg$enrichRes@result, geneList, keytype = "Entrez",
                            gene_cutoff = lfcCutoff, top = 0, bottom = 15)
      ggsave(filename=file.path(out.dir_sub,"RRA/EnrichedGeneView_go.neg.png"), p2,
             units = "in", width = 10, height = 7)
      grid.arrange(p2, ncol = 1)
    }
    if(!is.null(kegg.pos$enrichRes)){
      p3 = EnrichedGeneView(enrichment=kegg.pos$enrichRes@result, geneList, keytype = "Entrez",
                            gene_cutoff = lfcCutoff, top = 15, bottom = 0)
      ggsave(filename = file.path(out.dir_sub, "RRA/EnrichedGeneView_kegg.pos.png"), p3,
           units = "in", width = 10, height = 7)
      grid.arrange(p3, ncol = 1)
    }
    if(!is.null(go.pos$enrichRes)){
      p4 = EnrichedGeneView(enrichment=go.pos$enrichRes@result, geneList, keytype = "Entrez",
                            gene_cutoff = lfcCutoff, top = 15, bottom = 0)
      ggsave(filename=file.path(out.dir_sub, "RRA/EnrichedGeneView_go.pos.png"), p4,
           units = "in", width = 10, height = 7)
      grid.arrange(p4, ncol = 1)
    }
  }
  dev.off()
}
