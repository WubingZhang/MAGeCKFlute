#' Enrichment analysis for selected treatment related genes
#'
#' Do enrichment analysis for selected treatment related genes in 9-squares
#'
#' @docType methods
#' @name EnrichSquare
#' @rdname EnrichSquare
#'
#' @param beta Data frame, with columns of "Gene", "group", and "Diff".
#' @param id A character, indicating the gene column in the data.
#' @param keytype A character, "Symbol" or "Entrez".
#' @param x A character, indicating the x-axis in the 9-square scatter plot.
#' @param y A character, indicating the y-axis in the 9-square scatter plot.
#' @param enrich_method One of "ORT"(Over-Representing Test) and "HGT"(HyperGemetric test).
#' @param top An integer, specifying the number of pathways to show.
#' @param organism "hsa" or "mmu".
#' @param limit A two-length vector, specifying the min and
#' max size of pathways for enrichent analysis.
#' @param filename Suffix of output file name. NULL(default) means no output.
#' @param out.dir Path to save plot to (combined with filename).
#' @param width As in ggsave.
#' @param height As in ggsave.
#' @param verbose Boolean.
#' @param ... Other available parameters in ggsave.
#'
#' @return A list containing enrichment results for each group genes.
#' Each item in the returned list has two sub items:
#' \item{gridPlot}{an object created by \code{ggplot}, which can be assigned and further customized.}
#' \item{enrichRes}{a enrichResult instance.}
#'
#' @author Wubing Zhang

EnrichSquare <- function(beta, id = "Gene", keytype = "Symbol",
                         x = "Control", y = "Treatment",
                         enrich_method = "HGT", top = 5,
                         organism = "hsa", limit = c(2, 100),
                         filename=NULL, out.dir=".",
                         width=6.5, height=4, verbose = TRUE, ...){

  if(verbose) message(Sys.time(), " # Enrichment analysis of 9 Square grouped genes ...")
  ## ===========Enrichment===================
  gg = as.data.frame(beta[, c(id, x, y, "group")], stringsAsFactors = FALSE)
  colnames(gg) = c("Gene", "Control", "Treatment", "group")
  gg = gg[!(is.na(gg$Gene)|duplicated(gg$Gene)), ]
  gg$Diff = gg$Treatment - gg$Control
  idx1 = gg$group=="midleft"
  idx2 = gg$group=="topcenter"
  idx3 = gg$group=="midright"
  idx4 = gg$group=="bottomcenter"
  universe = gg$Gene

  #====GO_KEGG_enrichment=====
  geneList = gg$Control[idx1]; names(geneList) = gg$Gene[idx1]
  enrich1 = EnrichAnalyzer(geneList = geneList, universe=universe,
                           method = enrich_method,
                           type = "KEGG+REACTOME+GOBP+Complex",
                           organism=organism, pvalueCutoff = 1,
                           limit = limit, keytype = keytype, verbose = verbose)
  if(!is.null(enrich1) && nrow(enrich1@result)>0){
    kegg1 = enrich1@result[grepl("KEGG", enrich1@result$ID), ]
    gobp1 = enrich1@result[grepl("^GO", enrich1@result$ID), ]
    reactome1 = enrich1@result[grepl("REACTOME", enrich1@result$ID), ]
    complex1 = enrich1@result[grepl("CPX|CORUM", enrich1@result$ID), ]
    kegg1 = list(enrichRes = kegg1, gridPlot = EnrichedView(kegg1, top = 0, bottom = top)
                 + labs(title = "KEGG: Group1"))
    gobp1 = list(enrichRes = gobp1, gridPlot = EnrichedView(gobp1, top = 0, bottom = top)
                 + labs(title = "GOBP: Group1"))
    reactome1 = list(enrichRes = reactome1, gridPlot = EnrichedView(reactome1, top = 0, bottom = top)
                     + labs(title = "REACTOME: Group1"))
    complex1 = list(enrichRes = complex1, gridPlot = EnrichedView(complex1, top = 0, bottom = top)
                    + labs(title = "Complex: Group1"))
  }else{
    kegg1 = gobp1 = reactome1 = complex1 = list(enrichRes = NULL, gridPlot = noEnrichPlot())
  }

  #====GO_KEGG_enrichment=====
  geneList = gg$Treatment[idx2]; names(geneList) = gg$Gene[idx2]
  enrich2 = EnrichAnalyzer(geneList = geneList, universe=universe,
                           method = enrich_method, type = "KEGG+REACTOME+GOBP+Complex",
                           organism=organism, pvalueCutoff = 1,
                           limit = limit, keytype = keytype, verbose = verbose)
  if(!is.null(enrich2) && nrow(enrich2@result)>0){
    kegg2 = enrich2@result[grepl("KEGG", enrich2@result$ID), ]
    gobp2 = enrich2@result[grepl("^GO", enrich2@result$ID), ]
    reactome2 = enrich2@result[grepl("REACTOME", enrich2@result$ID), ]
    complex2 = enrich2@result[grepl("CPX|CORUM", enrich2@result$ID), ]
    kegg2 = list(enrichRes = kegg2, gridPlot = EnrichedView(kegg2, top = top, bottom = 0)
                 + labs(title = "KEGG: Group2"))
    gobp2 = list(enrichRes = gobp2, gridPlot = EnrichedView(gobp2, top = top, bottom = 0)
                 + labs(title = "GOBP: Group2"))
    reactome2 = list(enrichRes = reactome2, gridPlot = EnrichedView(reactome2, top = top, bottom = 0)
                     + labs(title = "REACTOME: Group2"))
    complex2 = list(enrichRes = complex2, gridPlot = EnrichedView(complex2, top = top, bottom = 0)
                    + labs(title = "Complex: Group2"))
  }else{
    kegg2 = gobp2 = reactome2 = complex2 = list(enrichRes = NULL, gridPlot = noEnrichPlot())
  }

  #====GO_KEGG_enrichment=====
  geneList = gg$Control[idx3]; names(geneList) = gg$Gene[idx3]
  enrich3 = EnrichAnalyzer(geneList = geneList, universe=universe,
                           method = enrich_method, type = "KEGG+REACTOME+GOBP+Complex",
                           organism=organism, pvalueCutoff = 1,
                           limit = limit, keytype = keytype, verbose = verbose)
  if(!is.null(enrich3) && nrow(enrich3@result)>0){
    kegg3 = enrich3@result[grepl("KEGG", enrich3@result$ID), ]
    gobp3 = enrich3@result[grepl("^GO", enrich3@result$ID), ]
    reactome3 = enrich3@result[grepl("REACTOME", enrich3@result$ID), ]
    complex3 = enrich3@result[grepl("CPX|CORUM", enrich3@result$ID), ]
    kegg3 = list(enrichRes = kegg3, gridPlot = EnrichedView(kegg3, top = top, bottom = 0)
                 + labs(title = "KEGG: Group3"))
    gobp3 = list(enrichRes = gobp3, gridPlot = EnrichedView(gobp3, top = top, bottom = 0)
                 + labs(title = "GOBP: Group3"))
    reactome3 = list(enrichRes = reactome3, gridPlot = EnrichedView(reactome3, top = top, bottom = 0)
                     + labs(title = "REACTOME: Group3"))
    complex3 = list(enrichRes = complex3, gridPlot = EnrichedView(complex3, top = top, bottom = 0)
                    + labs(title = "Complex: Group3"))
  }else{
    kegg3 = gobp3 = reactome3 = complex3 = list(enrichRes = NULL, gridPlot = noEnrichPlot())
  }

  #====GO_KEGG_enrichment=====
  geneList = gg$Treatment[idx4]; names(geneList) = gg$Gene[idx4]
  enrich4 = EnrichAnalyzer(geneList = geneList, universe=universe,
                           method = enrich_method, type = "KEGG+REACTOME+GOBP+Complex",
                           organism=organism, pvalueCutoff = 1,
                           limit = limit, keytype = keytype, verbose = verbose)
  if(!is.null(enrich4) && nrow(enrich4@result)>0){
    kegg4 = enrich4@result[grepl("KEGG", enrich4@result$ID), ]
    gobp4 = enrich4@result[grepl("^GO", enrich4@result$ID), ]
    reactome4 = enrich4@result[grepl("REACTOME", enrich4@result$ID), ]
    complex4 = enrich4@result[grepl("CPX|CORUM", enrich4@result$ID), ]
    kegg4 = list(enrichRes = kegg4, gridPlot = EnrichedView(kegg4, top = 0, bottom = top)
                 + labs(title = "KEGG: Group4"))
    gobp4 = list(enrichRes = gobp4, gridPlot = EnrichedView(gobp4, top = 0, bottom = top)
                 + labs(title = "GOBP: Group4"))
    reactome4 = list(enrichRes = reactome4, gridPlot = EnrichedView(reactome4, top = 0, bottom = top)
                     + labs(title = "REACTOME: Group4"))
    complex4 = list(enrichRes = complex4, gridPlot = EnrichedView(complex4, top = 0, bottom = top)
                    + labs(title = "Complex: Group4"))
  }else{
    kegg4 = gobp4 = reactome4 = complex4 = list(enrichRes = NULL, gridPlot = noEnrichPlot())
  }

  #====GO_KEGG_enrichment=====
  geneList = abs(gg$Diff[idx1|idx2]); names(geneList) = gg$Gene[idx1|idx2]
  enrich12 = EnrichAnalyzer(geneList = geneList, universe=universe,
                           method = enrich_method, type = "KEGG+REACTOME+GOBP+Complex",
                           organism=organism, pvalueCutoff = 1,
                           limit = limit, keytype = keytype, verbose = verbose)
  if(!is.null(enrich12) && nrow(enrich12@result)>0){
    kegg12 = enrich12@result[grepl("KEGG", enrich12@result$ID), ]
    gobp12 = enrich12@result[grepl("^GO", enrich12@result$ID), ]
    reactome12 = enrich12@result[grepl("REACTOME", enrich12@result$ID), ]
    complex12 = enrich12@result[grepl("CPX|CORUM", enrich12@result$ID), ]
    kegg12 = list(enrichRes = kegg12, gridPlot = EnrichedView(kegg12, top = top, bottom = 0)
                 + labs(title = "KEGG: Group12"))
    gobp12 = list(enrichRes = gobp12, gridPlot = EnrichedView(gobp12, top = top, bottom = 0)
                 + labs(title = "GOBP: Group12"))
    reactome12 = list(enrichRes = reactome12, gridPlot = EnrichedView(reactome12, top = top, bottom = 0)
                     + labs(title = "REACTOME: Group12"))
    complex12 = list(enrichRes = complex12, gridPlot = EnrichedView(complex12, top = top, bottom = 0)
                    + labs(title = "Complex: Group12"))
  }else{
    kegg12 = gobp12 = reactome12 = complex12 = list(enrichRes = NULL, gridPlot = noEnrichPlot())
  }

  #====GO_KEGG_enrichment=====
  geneList = gg$Control[idx1|idx3]; names(geneList) = gg$Gene[idx1|idx3]
  enrich13 = EnrichAnalyzer(geneList = geneList, universe=universe,
                           method = enrich_method, type = "KEGG+REACTOME+GOBP+Complex",
                           organism=organism, pvalueCutoff = 1,
                           limit = limit, keytype = keytype, verbose = verbose)
  if(!is.null(enrich13) && nrow(enrich13@result)>0){
    kegg13 = enrich13@result[grepl("KEGG", enrich13@result$ID), ]
    gobp13 = enrich13@result[grepl("^GO", enrich13@result$ID), ]
    reactome13 = enrich13@result[grepl("REACTOME", enrich13@result$ID), ]
    complex13 = enrich13@result[grepl("CPX|CORUM", enrich13@result$ID), ]
    kegg13 = list(enrichRes = kegg13, gridPlot = EnrichedView(kegg13, top = top, bottom = top)
                 + labs(title = "KEGG: Group13"))
    gobp13 = list(enrichRes = gobp13, gridPlot = EnrichedView(gobp13, top = top, bottom = top)
                 + labs(title = "GOBP: Group13"))
    reactome13 = list(enrichRes = reactome13, gridPlot = EnrichedView(reactome13, top = top, bottom = top)
                     + labs(title = "REACTOME: Group13"))
    complex13 = list(enrichRes = complex13, gridPlot = EnrichedView(complex13, top = top, bottom = top)
                    + labs(title = "Complex: Group13"))
  }else{
    kegg13 = gobp13 = reactome13 = complex13 = list(enrichRes = NULL, gridPlot = noEnrichPlot())
  }

  #====GO_KEGG_enrichment=====
  geneList = gg$Treatment[idx4|idx2]; names(geneList) = gg$Gene[idx4|idx2]
  enrich24 = EnrichAnalyzer(geneList = geneList, universe=universe,
                           method = enrich_method, type = "KEGG+REACTOME+GOBP+Complex",
                           organism=organism, pvalueCutoff = 1,
                           limit = limit, keytype = keytype, verbose = verbose)
  if(!is.null(enrich24) && nrow(enrich24@result)>0){
    kegg24 = enrich24@result[grepl("KEGG", enrich24@result$ID), ]
    gobp24 = enrich24@result[grepl("^GO", enrich24@result$ID), ]
    reactome24 = enrich24@result[grepl("REACTOME", enrich24@result$ID), ]
    complex24 = enrich24@result[grepl("CPX|CORUM", enrich24@result$ID), ]
    kegg24 = list(enrichRes = kegg24, gridPlot = EnrichedView(kegg24, top = top, bottom = top)
                 + labs(title = "KEGG: Group24"))
    gobp24 = list(enrichRes = gobp24, gridPlot = EnrichedView(gobp24, top = top, bottom = top)
                 + labs(title = "GOBP: Group24"))
    reactome24 = list(enrichRes = reactome24, gridPlot = EnrichedView(reactome24, top = top, bottom = top)
                     + labs(title = "REACTOME: Group24"))
    complex24 = list(enrichRes = complex24, gridPlot = EnrichedView(complex24, top = top, bottom = top)
                    + labs(title = "Complex: Group24"))
  }else{
    kegg24 = gobp24 = reactome24 = complex24 = list(enrichRes = NULL, gridPlot = noEnrichPlot())
  }

  #====GO_KEGG_enrichment=====
  geneList = gg$Diff[idx3|idx4]; names(geneList) = gg$Gene[idx3|idx4]
  enrich34 = EnrichAnalyzer(geneList = geneList, universe=universe,
                           method = enrich_method, type = "KEGG+REACTOME+GOBP+Complex",
                           organism=organism, pvalueCutoff = 1,
                           limit = limit, keytype = keytype, verbose = verbose)
  if(!is.null(enrich34) && nrow(enrich34@result)>0){
    kegg34 = enrich34@result[grepl("KEGG", enrich34@result$ID), ]
    gobp34 = enrich34@result[grepl("^GO", enrich34@result$ID), ]
    reactome34 = enrich34@result[grepl("REACTOME", enrich34@result$ID), ]
    complex34 = enrich34@result[grepl("CPX|CORUM", enrich34@result$ID), ]
    kegg34 = list(enrichRes = kegg34, gridPlot = EnrichedView(kegg34, top = 0, bottom = top)
                 + labs(title = "KEGG: Group34"))
    gobp34 = list(enrichRes = gobp34, gridPlot = EnrichedView(gobp34, top = 0, bottom = top)
                 + labs(title = "GOBP: Group34"))
    reactome34 = list(enrichRes = reactome34, gridPlot = EnrichedView(reactome34, top = 0, bottom = top)
                     + labs(title = "REACTOME: Group34"))
    complex34 = list(enrichRes = complex34, gridPlot = EnrichedView(complex34, top = 0, bottom = top)
                    + labs(title = "Complex: Group34"))
  }else{
    kegg34 = gobp34 = reactome34 = complex34 = list(enrichRes = NULL, gridPlot = noEnrichPlot())
  }

  ### Output results
  if(!is.null(filename)){
    if(!is.null(enrich1) && nrow(enrich1@result)>0){
      write.table(kegg1$enrichRes, file.path(out.dir,paste0("Group1_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(kegg1$gridPlot, filename=file.path(out.dir,paste0("Group1_kegg_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(reactome1$enrichRes, file.path(out.dir,paste0("Group1_reactome_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(reactome1$gridPlot, filename=file.path(out.dir,paste0("Group1_reactome_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(gobp1$enrichRes, file.path(out.dir,paste0("Group1_gobp_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(gobp1$gridPlot, filename=file.path(out.dir,paste0("Group1_gobp_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(complex1$enrichRes, file.path(out.dir,paste0("Group1_complex_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(complex1$gridPlot, filename=file.path(out.dir,paste0("Group1_complex_", filename,".png")),
             units = "in", width=6.5, height=4)
    }
    if(!is.null(enrich2) && nrow(enrich2@result)>0){
      write.table(kegg2$enrichRes, file.path(out.dir,paste0("Group2_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(kegg2$gridPlot, filename=file.path(out.dir,paste0("Group2_kegg_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(reactome2$enrichRes, file.path(out.dir,paste0("Group2_reactome_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(reactome2$gridPlot, filename=file.path(out.dir,paste0("Group2_reactome_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(gobp2$enrichRes, file.path(out.dir,paste0("Group2_gobp_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(gobp2$gridPlot, filename=file.path(out.dir,paste0("Group2_gobp_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(complex2$enrichRes, file.path(out.dir,paste0("Group2_complex_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(complex2$gridPlot, filename=file.path(out.dir,paste0("Group2_complex_", filename,".png")),
             units = "in", width=6.5, height=4)
    }
    if(!is.null(enrich3) && nrow(enrich3@result)>0){
      write.table(kegg3$enrichRes, file.path(out.dir,paste0("Group3_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(kegg3$gridPlot, filename=file.path(out.dir,paste0("Group3_kegg_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(reactome3$enrichRes, file.path(out.dir,paste0("Group3_reactome_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(reactome3$gridPlot, filename=file.path(out.dir,paste0("Group3_reactome_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(gobp3$enrichRes, file.path(out.dir,paste0("Group3_gobp_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(gobp3$gridPlot, filename=file.path(out.dir,paste0("Group3_gobp_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(complex3$enrichRes, file.path(out.dir,paste0("Group3_complex_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(complex3$gridPlot, filename=file.path(out.dir,paste0("Group3_complex_", filename,".png")),
             units = "in", width=6.5, height=4)
    }
    if(!is.null(enrich4) && nrow(enrich4@result)>0){
      write.table(kegg4$enrichRes, file.path(out.dir,paste0("Group4_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(kegg4$gridPlot, filename=file.path(out.dir,paste0("Group4_kegg_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(reactome4$enrichRes, file.path(out.dir,paste0("Group4_reactome_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(reactome4$gridPlot, filename=file.path(out.dir,paste0("Group4_reactome_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(gobp4$enrichRes, file.path(out.dir,paste0("Group4_gobp_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(gobp4$gridPlot, filename=file.path(out.dir,paste0("Group4_gobp_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(complex4$enrichRes, file.path(out.dir,paste0("Group4_complex_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(complex4$gridPlot, filename=file.path(out.dir,paste0("Group4_complex_", filename,".png")),
             units = "in", width=6.5, height=4)
    }
    if(!is.null(enrich13) && nrow(enrich13@result)>0){
      write.table(kegg13$enrichRes, file.path(out.dir,paste0("Group13_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(kegg13$gridPlot, filename=file.path(out.dir,paste0("Group13_kegg_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(reactome13$enrichRes, file.path(out.dir,paste0("Group13_reactome_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(reactome13$gridPlot, filename=file.path(out.dir,paste0("Group13_reactome_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(gobp13$enrichRes, file.path(out.dir,paste0("Group13_gobp_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(gobp13$gridPlot, filename=file.path(out.dir,paste0("Group13_gobp_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(complex13$enrichRes, file.path(out.dir,paste0("Group13_complex_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(complex13$gridPlot, filename=file.path(out.dir,paste0("Group13_complex_", filename,".png")),
             units = "in", width=6.5, height=4)
    }
    if(!is.null(enrich12) && nrow(enrich12@result)>0){
      write.table(kegg12$enrichRes, file.path(out.dir,paste0("Group12_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(kegg12$gridPlot, filename=file.path(out.dir,paste0("Group12_kegg_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(reactome12$enrichRes, file.path(out.dir,paste0("Group12_reactome_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(reactome12$gridPlot, filename=file.path(out.dir,paste0("Group12_reactome_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(gobp12$enrichRes, file.path(out.dir,paste0("Group12_gobp_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(gobp12$gridPlot, filename=file.path(out.dir,paste0("Group12_gobp_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(complex12$enrichRes, file.path(out.dir,paste0("Group12_complex_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(complex12$gridPlot, filename=file.path(out.dir,paste0("Group12_complex_", filename,".png")),
             units = "in", width=6.5, height=4)
    }
    if(!is.null(enrich24) && nrow(enrich24@result)>0){
      write.table(kegg24$enrichRes, file.path(out.dir,paste0("Group24_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(kegg24$gridPlot, filename=file.path(out.dir,paste0("Group24_kegg_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(reactome24$enrichRes, file.path(out.dir,paste0("Group24_reactome_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(reactome24$gridPlot, filename=file.path(out.dir,paste0("Group24_reactome_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(gobp24$enrichRes, file.path(out.dir,paste0("Group24_gobp_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(gobp24$gridPlot, filename=file.path(out.dir,paste0("Group24_gobp_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(complex24$enrichRes, file.path(out.dir,paste0("Group24_complex_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(complex24$gridPlot, filename=file.path(out.dir,paste0("Group24_complex_", filename,".png")),
             units = "in", width=6.5, height=4)
    }
    if(!is.null(enrich34) && nrow(enrich34@result)>0){
      write.table(kegg34$enrichRes, file.path(out.dir,paste0("Group34_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(kegg34$gridPlot, filename=file.path(out.dir,paste0("Group34_kegg_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(reactome34$enrichRes, file.path(out.dir,paste0("Group34_reactome_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(reactome34$gridPlot, filename=file.path(out.dir,paste0("Group34_reactome_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(gobp34$enrichRes, file.path(out.dir,paste0("Group34_gobp_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(gobp34$gridPlot, filename=file.path(out.dir,paste0("Group34_gobp_", filename,".png")),
             units = "in", width=6.5, height=4)
      write.table(complex34$enrichRes, file.path(out.dir,paste0("Group34_complex_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(complex34$gridPlot, filename=file.path(out.dir,paste0("Group34_complex_", filename,".png")),
             units = "in", width=6.5, height=4)
    }
  }
  return(list(kegg1=kegg1, gobp1=gobp1, reactome1=reactome1, complex1=complex1,
              kegg2=kegg2, gobp2=gobp2, reactome2=reactome2, complex2=complex2,
              kegg3=kegg3, gobp3=gobp3, reactome3=reactome3, complex3=complex3,
              kegg4=kegg4, gobp4=gobp4, reactome4=reactome4, complex4=complex4,
              kegg12=kegg12, gobp12=gobp12, reactome12=reactome12, complex12=complex12,
              kegg13=kegg13, gobp13=gobp13, reactome13=reactome13, complex13=complex13,
              kegg24=kegg24, gobp24=gobp24, reactome24=reactome24, complex24=complex24,
              kegg34=kegg34, gobp34=gobp34, reactome34=reactome34, complex34=complex34))
}

