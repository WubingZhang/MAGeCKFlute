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
#' @param score "lfc" (default) or "rra", specifying the score type.
#' @param cutoff A two-length vector (default: c(-1, 1)), specifying the cutoff
#' for negative selection and positive selection.
#' @param organism "hsa" or "mmu".
#' @param limit A two-length vector, specifying the minimal and
#' maximal size of gene sets for enrichent analysis.
#' @param pvalueCutoff A numeric, specifying pvalue cutoff of enrichment analysis, default 1.
#' @param proj A character, indicating the prefix of output file name.
#' @param width The width of summary pdf in inches.
#' @param height The height of summary pdf in inches.
#' @param outdir Output directory on disk.
#'
#' @author Wubing Zhang
#'
#' @return  All of the pipeline results is output into the \code{out.dir}/\code{proj}_Results,
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
#'     FluteRRA(rra.gene_summary, rra.sgrna_summary, proj="RRA", organism="hsa")
#' }
#'
#'
#' @export

FluteRRA <- function(gene_summary,
                     sgrna_summary = gsub("gene_summary", "sgrna_summary", gene_summary),
                     keytype = "Symbol",
                     score = c("lfc", "rra")[1],
                     cutoff = c(-1, 1),
                     organism = "hsa",
                     top = 5, toplabels = NULL,
                     limit = c(2, 200),
                     pvalueCutoff = 0.25,
                     proj = NA,
                     width = 12, height = 6,
                     outdir = "."){

  ## Prepare the output environment ##
  message(Sys.time(), " # Create output dir and pdf file ...")
  outdir = file.path(outdir, paste0("MAGeCKFlute_", proj))
  dir.create(file.path(outdir), showWarnings = FALSE)
  dir.create(file.path(outdir,"RRA"), showWarnings=FALSE)
  output_pdf = paste0("FluteRRA_", proj, "_", norm_method, ".pdf")

  pdf(file.path(outdir, output_pdf), width = width, height = height)

  ## Visualize the top essential genes ##
  message(Sys.time(), " # Read RRA result ...")
  dd = ReadRRA(gene_summary)
  dd.sgrna = ReadsgRRA(sgrna_summary)
  dd$LogFDR = -log10(dd$FDR)
  dd$Symbol = TransGeneID(dd$id, keytype, "Symbol", organism = organism)
  dd$EntrezID = TransGeneID(dd$id, keytype, "Entrez", organism = organism)
  dd.sgrna$Symbol = TransGeneID(dd.sgrna$Gene, keytype, "Entrez", organism = organism)
  idx1 = is.na(dd$EntrezID)
  idx2 = !is.na(dd$EntrezID) & duplicated(dd$EntrezID)
  idx = idx1|idx2
  if(sum(idx1)>0) warning(sum(idx1), " genes are not eligible: ",
                          paste0(dd$id[idx1], collapse = ", "))
  if(sum(idx2)>0) warning(sum(idx2), " genes have duplicate entrez id: ",
                          paste0(dd$id[idx2], collapse = ", "))
  dd = dd[!idx, ]

  p1 = ScatterView(dd, x = "Score", y = "LogFDR", label = "Symbol",
                   model = "volcano", top = top, toplabels = toplabels,
                   display_cut = TRUE)
  ggsave(file.path(outdir,"RRA/VolcanoView_RRA.png"), p1,
         units = "in", width = 6.5, height = 4)

  geneList = dd$Score
  names(geneList) = dd$Symbol
  dd$Rank = rank(dd$Score)
  p2 = ScatterView(dd, x = "Score", y = "Rank", label = "Symbol",
                   top = top, toplabels = toplabels, model = "rank",
                   display_cut = TRUE)
  ggsave(file.path(outdir,"RRA/RankView_Gene.png"), p2,
         units = "in", width = 6.5, height = 4)
  p3 = sgRankView(dd.sgrna, top = top, bottom = 0)
  ggsave(file.path(outdir,"RRA/RankView_sgRNA_top_.png"),
         p3, units = "in", width = 6.5, height = 4)
  p4 = sgRankView(dd.sgrna, top = 0, bottom = top)
  ggsave(file.path(outdir,"RRA/RankView_sgRNA_bottom_.png"),
         p4, units = "in", width = 6.5, height = 4)
  grid.arrange(p1, p2, p3, p4, ncol = 2)

  ## Enrichment analysis ##
  universe = dd$EntrezID
  geneList = dd$Score; names(geneList) = dd$EntrezID
  idx1 = dd$Score<cutoff[1]; idx2 = dd$Score>cutoff[2]
  kegg.neg = EnrichAnalyzer(geneList=geneList[idx1], universe=universe,
                            organism=organism, pvalueCutoff=pvalueCutoff,
                            limit = limit, keytype = "Entrez")
  p1 = EnrichedView(kegg.neg, bottom = top) + labs(title = "Enrichment: neg")

  kegg.pos = EnrichAnalyzer(geneList=geneList[idx2], universe=universe,
                            organism=organism, pvalueCutoff=pvalueCutoff,
                            limit = limit, keytype = "Entrez")
  p2 = EnrichedView(kegg.pos, top = top) + labs(title = "Enrichment: pos")
  grid.arrange(p1, p2, ncol = 1)

  ## Save enrichment results ##
  ggsave(file.path(outdir, "RRA/EnrichedView_neg.png"), p1,
         units = "in", width = 6.5, height = 4)
  ggsave(file.path(outdir, "RRA/EnrichedView_pos.png"), p2,
         units = "in", width = 6.5, height = 4)
  saveRDS(kegg.neg, file.path(outdir, "RRA/EnrichRes_kegg.neg.rds"))
  saveRDS(kegg.pos, file.path(outdir, "RRA/EnrichRes_kegg.pos.rds"))
  dev.off()
}
