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
#' @param keytype "Entrez" or "Symbol".
#' @param organism "hsa" or "mmu".
#' @param incorporateDepmap Boolean, indicating whether incorporate Depmap data into analysis.
#' @param cell_lines A character vector, specifying the cell lines in Depmap to be considered.
#' @param lineages A character vector, specifying the lineages in Depmap to be considered.
#' @param omitEssential Boolean, indicating whether omit common essential genes from the downstream analysis.
#' @param top An integer, specifying the number of top selected genes to be labeled
#' in rank figure and the number of top pathways to be shown.
#' @param toplabels A character vector, specifying interested genes to be labeled in rank figure.
#' @param scale_cutoff Boolean or numeric, specifying how many standard deviation will be used as cutoff.
#' @param limit A two-length vector, specifying the minimal and
#' maximal size of gene sets for enrichent analysis.
#' @param proj A character, indicating the prefix of output file name.
#' @param width The width of summary pdf in inches.
#' @param height The height of summary pdf in inches.
#' @param outdir Output directory on disk.
#' @param verbose Boolean
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
#' file1 = file.path(system.file("extdata", package = "MAGeCKFlute"),
#' "testdata/rra.gene_summary.txt")
#' file2 = file.path(system.file("extdata", package = "MAGeCKFlute"),
#'                   "testdata/rra.sgrna_summary.txt")
#' \dontrun{
#'     # Run the FluteRRA pipeline
#'     FluteRRA(file1, file2, proj="Pmel", organism="hsa", incorporateDepmap = FALSE,
#'     scale_cutoff = 1, outdir = "./")
#' }
#'
#' @export
#' @import ggplot2 gridExtra

FluteRRA <- function(gene_summary,
                     sgrna_summary = NULL,
                     keytype = "Symbol",
                     organism = "hsa",
                     incorporateDepmap = FALSE,
                     cell_lines = NA, lineages = "All",
                     omitEssential = TRUE,
                     top = 5, toplabels = NULL,
                     scale_cutoff = 2,
                     limit = c(2, 100),
                     proj = NA,
                     width = 12, height = 6,
                     outdir = ".",
                     verbose = TRUE){

  ## Prepare the output environment ##
  requireNamespace("ggplot2")
  message(Sys.time(), " # Create output dir and pdf file ...")
  outdir = file.path(outdir, paste0("MAGeCKFlute_", proj))
  dir.create(file.path(outdir), showWarnings = FALSE)
  dir.create(file.path(outdir,"RRA"), showWarnings=FALSE)
  output_pdf = paste0("FluteRRA_", proj, ".pdf")

  pdf(file.path(outdir, output_pdf), width = width, height = height)

  ## Visualize the top essential genes ##
  message(Sys.time(), " # Read RRA result ...")
  dd = ReadRRA(gene_summary)
  dd$LogFDR = -log10(dd$FDR)
  if(tolower(keytype)!="entrez"){
    dd$EntrezID = TransGeneID(dd$id, keytype, "Entrez", organism = organism)
  }else{
    dd$EntrezID = dd$id
  }
  if(tolower(keytype)!="symbol"){
    dd$Symbol = TransGeneID(dd$id, keytype, "Symbol", organism = organism)
  }else{
    dd$Symbol = dd$id
  }
  if(organism != "hsa"){
    dd$HumanGene = TransGeneID(dd$id, keytype, "Symbol",
                               fromOrg = organism, toOrg = "hsa")
  }else{
    dd$HumanGene = dd$Symbol
  }
  idx1 = is.na(dd$EntrezID)
  idx2 = !is.na(dd$EntrezID) & duplicated(dd$EntrezID)
  idx = idx1|idx2
  if(sum(idx1)>0) message(sum(idx1), " genes fail to convert into Entrez IDs: ",
                          paste0(dd$id[idx1], collapse = ", "))
  if(sum(idx2)>0) message(sum(idx2), " genes have duplicate Entrez IDs: ",
                          paste0(dd$id[idx2], collapse = ", "))
  dd = dd[!idx, ]

  if(!is.null(sgrna_summary)){
    dd.sgrna = ReadsgRRA(sgrna_summary)
    dd.sgrna = dd.sgrna[dd.sgrna$Gene%in%dd$id, ]
    if(tolower(keytype)!="entrez"){
      dd.sgrna$EntrezID = TransGeneID(dd.sgrna$id, keytype, "Entrez", organism = organism)
    }else{
      dd.sgrna$EntrezID = dd.sgrna$id
    }
    if(tolower(keytype)!="symbol"){
      dd.sgrna$Symbol = TransGeneID(dd.sgrna$id, keytype, "Symbol", organism = organism)
    }else{
      dd.sgrna$Symbol = dd.sgrna$id
    }
    if(organism != "hsa"){
      dd.sgrna$HumanGene = TransGeneID(dd.sgrna$Gene, keytype, "Symbol",
                                       fromOrg = organism, toOrg = "hsa")
    }else{
      dd.sgrna$HumanGene = dd.sgrna$Symbol
    }
  }else{
    dd.sgrna = data.frame(sgrna = NA, Gene = NA, LFC = NA,
                          EntrezID = NA, Symbol = NA, HumanGene = NA)
  }
  cutoff = c(-CutoffCalling(dd$Score, scale_cutoff),
             CutoffCalling(dd$Score, scale_cutoff))
  write.table(dd, file.path(outdir, paste0("RRA/", proj, "_processed_data.txt")),
              sep = "\t", row.names = FALSE, quote = FALSE)

  if(incorporateDepmap){
    dd = IncorporateDepmap(dd, symbol = "HumanGene", cell_lines = cell_lines,
                           lineages = lineages)
    ## Nine-squares ##
    p.square = ScatterView(dd, x = "Depmap", y = "Score", label = "Symbol",
                           groups = c("midleft", "topcenter", "midright", "bottomcenter"),
                           groupnames = c("Group1", "Group2", "Group3", "Group4"),
                           auto_cut_x = TRUE, y_cut = cutoff,
                           auto_cut_diag = TRUE, top = top,
                           display_cut = TRUE) + ylab("Customized")
    ggsave(file.path(outdir, "RRA/SquareView_Customized_Depmap.png"),
           p.square, width = 5, height = 4)
    write.table(dd, file.path(outdir, paste0("RRA/", proj, "_incorporate_depmap.txt")),
                sep = "\t", row.names = FALSE, quote = FALSE)

    E1 = EnrichSquare(p.square$data, id = "EntrezID", keytype = "entrez",
                      x = "Depmap", y = "Score", top = top,
                      organism=organism, filename="RRA", limit = limit,
                      out.dir=file.path(outdir, "RRA/"))
  }
  if(omitEssential){
    dd = OmitCommonEssential(dd, symbol = "HumanGene")
    dd.sgrna = OmitCommonEssential(dd.sgrna, symbol = "HumanGene")
    write.table(dd, file.path(outdir, paste0("RRA/", proj, "_omit_essential.txt")),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
  dd$RandomIndex = sample(1:nrow(dd), nrow(dd))
  p1 = ScatterView(dd, x = "Score", y = "LogFDR", label = "Symbol",
                   x_cut = cutoff, model = "volcano", top = top,
                   toplabels = toplabels, display_cut = TRUE)
  ggsave(file.path(outdir,"RRA/VolcanoView_RRA.png"), p1,
         units = "in", width = 5, height = 4)

  # geneList = dd$Score
  # names(geneList) = dd$Symbol
  # dd$Rank = rank(dd$Score)
  # p2 = ScatterView(dd, x = "Rank", y = "Score", label = "Symbol",
  #                  y_cut = cutoff, groups = c("top", "bottom"),
  #                  top = top, toplabels = toplabels,
  #                  display_cut = TRUE)
  # ggsave(file.path(outdir,"RRA/RankView_Gene.png"), p2,
  #        units = "in", width = 3, height = 5)
  p2 = sgRankView(dd.sgrna, top = top, bottom = top)
  ggsave(file.path(outdir,"RRA/RankView_sgRNA.png"), p2, units = "in", width = 6.5, height = 5)
  p3 = ScatterView(dd[dd$Score>0, ], x = "RandomIndex", y = "Score", label = "Symbol",
                   y_cut = cutoff, groups = "top", top = top)
  ggsave(file.path(outdir, "RRA/ScatterView_Positive.png"), p3, width = 5, height = 4)
  p4 = ScatterView(dd[dd$Score<0, ], x = "RandomIndex", y = "Score", label = "Symbol",
                   auto_cut_y = TRUE, groups = "bottom", top = top)
  ggsave(file.path(outdir, "RRA/ScatterView_Negative.png"), p4, width = 5, height = 4)
  gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)

  ## Enrichment analysis ##
  # dd = dd[!(is.na(dd$HumanGene)|duplicated(dd$HumanGene)), ]
  universe = dd$HumanGene
  geneList = dd$Score; names(geneList) = dd$HumanGene
  enrichRes = EnrichAnalyzer(geneList=geneList, universe=universe,
                             organism="hsa", pvalueCutoff=1,
                             method = "GSEA", limit = limit,
                             type = "KEGG+REACTOME+GOBP+Complex")
  kegg.pos = enrichRes@result[enrichRes@result$NES>0, ]
  if(!is.null(kegg.pos) && nrow(kegg.pos)>0){
    keggA = kegg.pos[grepl("KEGG", kegg.pos$ID), ]
    gobpA = kegg.pos[grepl("^GO", kegg.pos$ID), ]
    reactomeA = kegg.pos[grepl("REACTOME", kegg.pos$ID), ]
    complexA = kegg.pos[grepl("CPX|CORUM", kegg.pos$ID), ]
    keggA = list(enrichRes = keggA, gridPlot = EnrichedView(keggA, top = top, bottom = 0)
                 + labs(title = "KEGG: positive"))
    gobpA = list(enrichRes = gobpA, gridPlot = EnrichedView(gobpA, top = top, bottom = 0)
                 + labs(title = "GOBP: positive"))
    reactomeA = list(enrichRes = reactomeA, gridPlot = EnrichedView(reactomeA, top = top, bottom = 0)
                     + labs(title = "REACTOME: positive"))
    complexA = list(enrichRes = complexA, gridPlot = EnrichedView(complexA, top = top, bottom = 0)
                    + labs(title = "Complex: positive"))
  }else{
    keggA = gobpA = reactomeA = complexA = list(enrichRes = NULL, gridPlot = noEnrichPlot())
  }
  kegg.neg = enrichRes@result[enrichRes@result$NES<0, ]
  if(!is.null(kegg.neg) && nrow(kegg.neg)>0){
    keggB = kegg.neg[grepl("KEGG", kegg.neg$ID), ]
    gobpB = kegg.neg[grepl("^GO", kegg.neg$ID), ]
    reactomeB = kegg.neg[grepl("REACTOME", kegg.neg$ID), ]
    complexB = kegg.neg[grepl("CPX|CORUM", kegg.neg$ID), ]
    keggB = list(enrichRes = keggB,
                 gridPlot = EnrichedView(keggB, top = 0, bottom = top)
                 + labs(title = "KEGG: negative"))
    gobpB = list(enrichRes = gobpB,
                 gridPlot = EnrichedView(gobpB, top = 0, bottom = top)
                 + labs(title = "GOBP: negative"))
    reactomeB = list(enrichRes = reactomeB,
                     gridPlot = EnrichedView(reactomeB, top = 0, bottom = top)
                     + labs(title = "REACTOME: negative"))
    complexB = list(enrichRes = complexB,
                    gridPlot = EnrichedView(complexB, top = 0, bottom = top)
                    + labs(title = "Complex: negative"))
  }else{
    keggB = gobpB = reactomeB = complexB = list(enrichRes = NULL, gridPlot = noEnrichPlot())
  }
  gridExtra::grid.arrange(keggA$gridPlot, gobpA$gridPlot, reactomeA$gridPlot, complexA$gridPlot, ncol = 2)
  gridExtra::grid.arrange(keggB$gridPlot, gobpB$gridPlot, reactomeB$gridPlot, complexB$gridPlot, ncol = 2)

  ## Save enrichment results ##
  if(!is.null(kegg.pos) && nrow(kegg.pos)>0){
    write.table(keggA$enrichRes, file.path(outdir, "RRA/Positive_kegg.txt"),
                sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
    ggsave(keggA$gridPlot, filename=file.path(outdir, "RRA/Positive_kegg.png"),
           units = "in", width=6.5, height=4)
    write.table(reactomeA$enrichRes, file.path(outdir, "RRA/Positive_reactome.txt"),
                sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
    ggsave(reactomeA$gridPlot, filename=file.path(outdir, "RRA/Positive_reactome.png"),
           units = "in", width=6.5, height=4)
    write.table(gobpA$enrichRes, file.path(outdir, "RRA/Positive_gobp.txt"),
                sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
    ggsave(gobpA$gridPlot, filename=file.path(outdir, "RRA/Positive_gobp.png"),
           units = "in", width=6.5, height=4)
    write.table(complexA$enrichRes, file.path(outdir, "RRA/Positive_complex.txt"),
                sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
    ggsave(complexA$gridPlot, filename=file.path(outdir, "RRA/Positive_complex.png"),
           units = "in", width=6.5, height=4)
  }
  if(!is.null(kegg.neg) && nrow(kegg.neg)>0){
    write.table(keggB$enrichRes, file.path(outdir, "RRA/Negative_kegg.txt"),
                sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
    ggsave(keggB$gridPlot, filename=file.path(outdir, "RRA/Negative_kegg.png"),
           units = "in", width=6.5, height=4)
    write.table(reactomeB$enrichRes, file.path(outdir, "RRA/Negative_reactome.txt"),
                sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
    ggsave(reactomeB$gridPlot, filename=file.path(outdir, "RRA/Negative_reactome.png"),
           units = "in", width=6.5, height=4)
    write.table(gobpB$enrichRes, file.path(outdir, "RRA/Negative_gobp.txt"),
                sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
    ggsave(gobpB$gridPlot, filename=file.path(outdir, "RRA/Negative_gobp.png"),
           units = "in", width=6.5, height=4)
    write.table(complexB$enrichRes, file.path(outdir, "RRA/Negative_complex.txt"),
                sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
    ggsave(complexB$gridPlot, filename=file.path(outdir, "RRA/Negative_complex.png"),
           units = "in", width=6.5, height=4)
  }
  if(incorporateDepmap){
    gridExtra::grid.arrange(p.square, ncol = 1)
    gridExtra::grid.arrange(E1$kegg1$gridPlot, E1$reactome1$gridPlot,
                 E1$gobp1$gridPlot, E1$complex1$gridPlot, ncol = 2)
    gridExtra::grid.arrange(E1$kegg2$gridPlot, E1$reactome2$gridPlot,
                 E1$gobp2$gridPlot, E1$complex2$gridPlot, ncol = 2)
    gridExtra::grid.arrange(E1$kegg3$gridPlot, E1$reactome3$gridPlot,
                 E1$gobp3$gridPlot, E1$complex3$gridPlot, ncol = 2)
    gridExtra::grid.arrange(E1$kegg4$gridPlot, E1$reactome4$gridPlot,
                 E1$gobp4$gridPlot, E1$complex4$gridPlot, ncol = 2)
    gridExtra::grid.arrange(E1$kegg12$gridPlot, E1$reactome12$gridPlot,
                 E1$gobp12$gridPlot, E1$complex12$gridPlot, ncol = 2)
    gridExtra::grid.arrange(E1$kegg13$gridPlot, E1$reactome13$gridPlot,
                 E1$gobp13$gridPlot, E1$complex13$gridPlot, ncol = 2)
    gridExtra::grid.arrange(E1$kegg24$gridPlot, E1$reactome24$gridPlot,
                 E1$gobp24$gridPlot, E1$complex24$gridPlot, ncol = 2)
    gridExtra::grid.arrange(E1$kegg34$gridPlot, E1$reactome34$gridPlot,
                 E1$gobp34$gridPlot, E1$complex34$gridPlot, ncol = 2)
  }
  dev.off()
}
