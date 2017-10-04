enrich.GOstats <- function(gene, universe=NULL, type=c("KEGG", "BP", "MF", "CC"), organism = "hsa",
                           pvalueCutoff = 0.05, pAdjustMethod = "none"){
  #======
  DS = toupper(type[1])
  over.sum=data.frame()
  orgdb = c("org.Hs.eg.db", "org.Mm.eg.db")
  names(orgdb) = c("hsa", "mmu")
  #========
  if(DS == "KEGG"){
    loginfo(paste('Running KEGG for list of entrezIDs'))
    params <- new("KEGGHyperGParams", categoryName="KEGG",
                  geneIds=gene, universeGeneIds=universe,
                  annotation=orgdb[organism], pvalueCutoff=pvalueCutoff,testDirection="over")
    #loginfo('    Starting HyperG test')
    over <- hyperGTest(params)
    #loginfo('    HyperG test done')
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
    loginfo(paste('Running GO', DS, 'for list of entrezIDs'))
    params <- new("GOHyperGParams", annotation=orgdb[organism],
                  geneIds = gene, universeGeneIds = universe,
                  ontology=DS, pvalueCutoff=pvalueCutoff, testDirection="over")
    over <- hyperGTest(params)
    over.sum <- summary(over)
  }

  #==================
  result = list()
  boo=!is.null(over.sum) & (nrow(over.sum)>0)
  if(boo){
    glist <- geneIdsByCategory(over)
    glist <- sapply(glist, function(.ids) {
      .sym <- .ids
      # .sym <- mget(.ids, envir=org.Hs.egSYMBOL, ifnotfound=NA)
      .sym[is.na(.sym)] <- .ids[is.na(.sym)]
      paste(.sym, collapse="/")
    })
    over.sum$geneID <- glist[as.character(over.sum[,1])]
    geneID = strsplit(over.sum$geneID, "/")
    geneName = lapply(geneID, function(gid){
      SYMBOL = TransGeneID(gid, "ENTREZID", "SYMBOL", organism)
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
