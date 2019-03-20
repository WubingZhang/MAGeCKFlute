#' Downstream analysis based on MAGeCK-MLE result
#'
#' Identify essential genes from MAGeCK MLE results
#'
#' @docType methods
#' @name subFluteMLE
#' @rdname subFluteMLE
#'
#' @param dd A data frame.
#' @param ctrlname A character vector, specifying the names of control samples in \code{dd}.
#' @param treatname A character vector, specifying the names of treatment samples in \code{dd}.
#' @param organism "hsa" or "mmu".
#'
#' @param scale_cutoff Boolean or numeric, whether scale cutoff to whole genome level,
#' or how many standard deviation will be used as cutoff.
#' @param top An integer, specifying number of top selected genes to be labeled in rank figure.
#' @param bottom An integer, specifying number of bottom selected genes to be labeled in rank figure.
#' @param interestGenes A character vector, specifying interested genes to be labeled in rank figure.
#'
#' @param limit A two-length vector (default: c(3, 50)), specifying the minimal and
#' maximal size of gene sets for enrichent analysis.
#' @param pvalueCutoff A numeric, specifying pvalue cutoff of enrichment analysis, default 1.
#' @param adjust One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param enrich_kegg One of "ORT" and "HGT"(HyperGemetric test).
#' @param posControl A character vector, specifying a list of positive control gene symbols.
#'
#' @author Wubing Zhang
#' @import ggplot2 stats grDevices utils gridExtra grid

function(dd, ctrlname, treatname, scale_cutoff, top, bottom, interestGenes,
         main, outdir){
  ## Distribution of all genes ##
  idx_distr = c(ctrlname, treatname)
  P1 = ViolinView(dd[, idx_distr], main = main)
  P2 = DensityView(dd[, idx_distr], main = main)
  P3 = DensityDiffView(dd, ctrlname, treatname, main = main)
  P4 = MAView(dd, ctrlname, treatname, main = main)

  ## Distribution of positive control genes ##
  data(Zuber_Essential)
  if(is.null(posControl))
    idx = toupper(dd$Gene) %in% toupper(Zuber_Essential$GeneSymbol)
  else
    idx = which(dd$Gene %in% posControl)
  P5 = ViolinView(dd[idx, idx_distr], ylab = "Essential.B.S.", main = main)
  P6 = DensityView(dd[idx, idx_distr], xlab = "Essential.B.S.", main = main)
  P7 = CellCycleView(dd[, idx_distr], ctrlname, treatname, main = main)
  P8 = CellCycleView(dd[idx, idx_distr], ctrlname, treatname, main = main)

  ## Identify essential genes based on deviation between Treat.BS and Ctrl.BS ##
  dd$Control = rowMeans(dd[, ctrlname, drop = FALSE])
  dd$Treatment = rowMeans(dd[, treatname, drop = FALSE])
  dd.diff = dd$Treatment - dd$Control
  names(dd.diff) = dd$Gene
  P9 = ScatterView(dd, main = main, scale_cutoff = scale_cutoff)
  P10 = RankView(dd.diff, genelist = interestGenes, top = top, bottom = bottom, main = main,
                cutoff = c(-CutoffCalling(dd.diff, scale = scale_cutoff),
                           CutoffCalling(dd.diff, scale = scale_cutoff)))

  ## groupA and groupB ##
  E2 = EnrichAB(P1$data, pvalue = pvalueCutoff, enrich_method = enrich_kegg,
                organism = organism, adjust = adjust, gsea = gsea, limit = limit)

  ## 9-Square model ##
  P11 = SquareView(dd, label="Gene", main = main, groups = c("midleft", "midright"),
                   groupnames = c("Group1", "Group3"))
  P12 = SquareView(dd, label="Gene", main = main, groups = c("topcenter", "bottomcenter"),
                   groupnames = c("Group2", "Group4"))
  P13 = SquareView(dd, label="Gene", main = main)
  E9 = EnrichSquare(P13$data, pvalue = pvalueCutoff, adjust = adjust,
                    enrich_method = enrich_kegg, organism = organism, limit = limit)
}
