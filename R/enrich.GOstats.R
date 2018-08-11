#' Do enrichment analysis using GOstats
#'
#' Do enrichment analysis using GOstats method
#'
#' @docType methods
#' @name enrich.GOstats
#' @rdname enrich.GOstats
#' @aliases enrichGOstats
#'
#' @param gene A character vector, specifying the genelist to do enrichment analysis.
#' @param universe A character vector, specifying the backgound genelist, default is whole genome.
#' @param type Geneset category for testing, KEGG(default).
#' @param organism A character, specifying organism, such as "hsa" or "Human"(default), and "mmu" or "Mouse"
#' @param pvalueCutoff Pvalue cutoff.
#' @param pAdjustMethod One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#'
#' @return A enrichResult instance.
#'
#' @author Wubing Zhang
#'
#' @seealso \code{\link{enrich.HGT}}
#' @seealso \code{\link{enrich.DAVID}}
#' @seealso \code{\link{enrich.GSE}}
#' @seealso \code{\link{enrich.ORT}}
#' @seealso \code{\link{enrichment_analysis}}
#' @seealso \code{\link[DOSE]{enrichResult-class}}
#'
#' @examples
#' data(geneList, package = "DOSE")
#' genes <- names(geneList)[1:100]
#' enrichRes <- enrich.GOstats(genes, type="BP")
#' head(enrichRes@result)
#'
#' @import GOstats Category
#' @import DOSE
#'
#' @export


enrich.GOstats <- function(gene, universe=NULL, type=c("KEGG", "BP", "MF", "CC"), organism = "hsa",
                           pvalueCutoff = 0.25, pAdjustMethod = "BH"){
  requireNamespace("GOstats", quietly = TRUE)
  requireNamespace("Category", quietly = TRUE)
  gene = unique(as.character(gene))
  if(!is.null(universe)) universe = unique(as.character(universe))
  #======
  DS = toupper(type[1])
  over.sum=data.frame()
  orgdb = getOrg(organism)$pkg
  #========
  if(DS == "KEGG"){
    params <- new("KEGGHyperGParams", categoryName="KEGG", geneIds=gene, universeGeneIds=universe,
                  annotation=orgdb, pvalueCutoff=pvalueCutoff,testDirection="over")
    over <- hyperGTest(params)
    over.sum <- summary(over)
#
#     glist <- geneIdsByCategory(over)
#     glist <- sapply(glist, function(.ids) {
#       .sym <- .ids
#       # .sym <- mget(.ids, envir=org.Hs.egSYMBOL, ifnotfound=NA)
#       .sym[is.na(.sym)] <- .ids[is.na(.sym)]
#       paste(.sym, collapse="/")
#     })
#     over.sum$geneID <- glist[as.character(over.sum$KEGGID)]
#     # over.sum$Symbols <- glist[as.character(over.sum$KEGGID)]
  }

  if(DS %in% c("BP", "CC", "MF")){
    # gene ontology background
    params <- new("GOHyperGParams", annotation=orgdb,
                  geneIds = gene, universeGeneIds = universe,
                  ontology=DS, pvalueCutoff=pvalueCutoff, testDirection="over")
    over <- hyperGTest(params)
    over.sum <- summary(over)
  }

  #==================
  result = list()
  if(!is.null(over.sum) & (nrow(over.sum)>0)){
    glist <- geneIdsByCategory(over)
    glist <- sapply(glist, function(.ids) {
      .sym <- .ids
      # .sym <- mget(.ids, envir=org.Hs.egSYMBOL, ifnotfound=NA)
      .sym[is.na(.sym)] <- .ids[is.na(.sym)]
      paste(.sym, collapse="/")
    })
    over.sum$geneID <- glist[as.character(over.sum[,1])]
    geneID = strsplit(over.sum$geneID, "/")
    allsymbol = TransGeneID(gene, "Entrez", "Symbol", organism = organism)
    geneName = lapply(geneID, function(gid){
      SYMBOL = allsymbol[gid]
      paste(SYMBOL, collapse = "/")
    })
    over.sum$geneName <- unlist(geneName)

    result$ID = over.sum[,1]
    result$Description = over.sum$Term
    result$GeneRatio = paste0(over.sum$Count, "/", over.sum$Size)
    result$BgRatio = paste0(over.sum$Size, "/", length(universe))
    result$pvalue = over.sum$Pvalue
    result$p.adjust = p.adjust(result$pvalue, method=pAdjustMethod)
    result$geneID = over.sum$geneID
    result$geneName = over.sum$geneName
    result$Count = over.sum$Count
    result = as.data.frame(result)
    idx = which(result$pvalue<=pvalueCutoff & result$p.adjust<=pvalueCutoff)
    result = result[idx, ]
    result = result[order(result$p.adjust),]
  }

  new("enrichResult",
      result         = result,
      pvalueCutoff   = pvalueCutoff,
      pAdjustMethod  = pAdjustMethod,
      organism       = organism,
      ontology       = type, ## as.character(x$Category[1]),
      gene           = as.character(gene),
      keytype        = "ENTREZID")
}
