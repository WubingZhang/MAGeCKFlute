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
#' @param limit A two-length vector, specifying the minimal and
#' maximal size of gene sets for enrichent analysis.
#' @param pvalueCutoff A numeric, specifying pvalue cutoff of enrichment analysis, default 1.
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
#' data("rra.gene_summary")
#' data("rra.sgrna_summary")
#' \dontrun{
#'     # Run the FluteRRA pipeline
#'     FluteRRA(rra.gene_summary, rra.sgrna_summary, prefix="RRA", organism="hsa")
#' }
#'
#'
#' @export

FluteRRA <- function(gene_summary, sgrna_summary,
                     lfcCutoff = c(-1, 1),
                     organism = "hsa",
                     limit = c(2, 200),
                     pvalueCutoff = 0.25,
                     prefix = "Test",
                     width = 12, height = 6,
                     outdir = "."){

  ## Prepare the output environment ##
  message(Sys.time(), " # Create output dir and pdf file ...")
  prefix = gsub(".*\\/|.*\\\\", "", prefix)
  out.dir_sub=file.path(outdir, paste0(prefix, "_Flute_Results"))
  dir.create(file.path(out.dir_sub), showWarnings=FALSE)
  dir.create(file.path(out.dir_sub,"RRA"), showWarnings=FALSE)

  output_pdf = paste0(prefix,"_Flute.rra_summary.pdf")
  pdf(file.path(out.dir_sub, output_pdf), width = width, height = height)

  ## Visualize the top essential genes ##
  message(Sys.time(), " # Read RRA result ...")
  dd = ReadRRA(gene_summary)
  dd.sgrna = ReadsgRRA(sgrna_summary)

  p1 = VolcanoView(dd, x = "LFC", y = "FDR", Label = "Official")
  ggsave(file.path(out.dir_sub,"RRA/VolcanoView_RRA.png"), p1,
         units = "in", width = 6.5, height = 4)

  geneList = dd$LFC
  names(geneList) = dd$Official
  p2 = RankView(geneList)
  p2 = p2 + labs(x = "Log2 Fold Change")
  ggsave(file.path(out.dir_sub,"RRA/RankView_Gene.png"), p2,
         units = "in", width = 6.5, height = 4)

  p3 = sgRankView(dd.sgrna, top = 8, bottom = 0)
  ggsave(file.path(out.dir_sub,"RRA/RankView_sgRNA_top8_.png"),
         p3, units = "in", width = 6.5, height = 4)
  p4 = sgRankView(dd.sgrna, top = 0, bottom = 8)
  ggsave(file.path(out.dir_sub,"RRA/RankView_sgRNA_bott8_.png"),
         p4, units = "in", width = 6.5, height = 4)
  grid.arrange(p1, p2, p3, p4, ncol = 2)

  ## Enrichment analysis ##
  dd$EntrezID = TransGeneID(dd$Official, "Symbol", "Entrez", organism = organism)
  idx = is.na(dd$EntrezID) | duplicated(dd$EntrezID)
  dd = dd[!idx,]

  universe = dd$EntrezID
  geneList = dd$LFC; names(geneList) = dd$EntrezID
  idx1 = dd$LFC<lfcCutoff[1]; idx2 = dd$LFC>lfcCutoff[2]
  kegg.neg = EnrichAnalyzer(geneList=geneList[idx1], universe=universe,
                            organism=organism, pvalueCutoff=pvalueCutoff,
                            limit = limit)
  p1 = EnrichedView(kegg.neg, bottom = 8) + labs(title = "Enrichment: neg")
  p2 = EnrichedGeneView(kegg.neg, geneList, keytype = "Entrez",
                        gene_cutoff = lfcCutoff, top = 0, bottom = 5)

  kegg.pos = EnrichAnalyzer(geneList=geneList[idx2], universe=universe,
                            organism=organism, pvalueCutoff=pvalueCutoff,
                            limit = limit)
  p3 = EnrichedView(kegg.pos, top = 8) + labs(title = "Enrichment: pos")
  p4 = EnrichedGeneView(kegg.pos, geneList, keytype = "Entrez",
                        gene_cutoff = lfcCutoff, top = 5, bottom = 0)
  grid.arrange(p2, p4, ncol = 1)

  ## Save enrichment results ##
  ggsave(file.path(out.dir_sub, "RRA/EnrichedView_neg.png"), p1,
         units = "in", width = 6.5, height = 4)
  ggsave(file.path(out.dir_sub,"RRA/EnrichedGeneView_neg.png"), p2,
         units = "in", width = 10, height = 7)
  ggsave(file.path(out.dir_sub, "RRA/EnrichedView_pos.png"), p3,
         units = "in", width = 6.5, height = 4)
  ggsave(file.path(out.dir_sub,"RRA/EnrichedGeneView_pos.png"), p4,
         units = "in", width = 10, height = 7)
  saveRDS(kegg.neg, file.path(out.dir_sub, "RRA/EnrichRes_kegg.neg.rds"))
  saveRDS(kegg.pos, file.path(out.dir_sub, "RRA/EnrichRes_kegg.pos.rds"))
  dev.off()
}
