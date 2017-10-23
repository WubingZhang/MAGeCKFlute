enrich.ORT <- function(gene, universe=NULL, type="KEGG", organism = "hsa",pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2, pAdjustMethod = "BH",minGSSize = 2, maxGSSize = 500){
  loginfo(paste('Running Over-Representation Test for list of entrezIDs'))
  gene = unique(as.character(gene))
  universe = unique(as.character(universe))
  orgdb = c("org.Hs.eg.db", "org.Mm.eg.db")
  names(orgdb) = c("hsa", "mmu")
  #=======================
  if(type %in% c("BP", "CC", "MF")){
    enrichedRes = enrichGO(gene=gene,  universe=universe,  ont = type, OrgDb=orgdb[organism],
                           pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff,
                           minGSSize=minGSSize, maxGSSize=maxGSSize)
    if(!is.null(enrichedRes))
        enrichedRes = simplify(enrichedRes, cutoff=0.7, by="p.adjust", select_fun=min)
  }
  if(type == "KEGG"){
    enrichedRes = enrichKEGG(gene=gene,  universe=universe, organism = organism,
                             pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff, qvalueCutoff=qvalueCutoff,
                             minGSSize = minGSSize, maxGSSize = maxGSSize, use_internal_data=FALSE)
  }
  if(type == "DO"){
    enrichedRes = enrichDO(gene=gene,  universe=universe, pAdjustMethod = pAdjustMethod,pvalueCutoff = pvalueCutoff,
                           qvalueCutoff=qvalueCutoff, minGSSize = minGSSize, maxGSSize = maxGSSize)
  }
  if(type == "MKEGG"){
    enrichedRes = enrichMKEGG(gene=gene,  universe=universe, organism = organism,pAdjustMethod=pAdjustMethod,
                              pvalueCutoff = pvalueCutoff,qvalueCutoff=qvalueCutoff, minGSSize = minGSSize, maxGSSize = maxGSSize)
  }
  if(type == "NCG"){
    enrichedRes = enrichNCG(gene=gene,  universe=universe, pAdjustMethod = pAdjustMethod,pvalueCutoff=pvalueCutoff,
                            qvalueCutoff=qvalueCutoff, minGSSize = minGSSize, maxGSSize = maxGSSize)
  }
  if(!is.null(enrichedRes)){
    geneID = strsplit(enrichedRes@result$geneID, "/")
    geneName = lapply(geneID, function(gid){
      SYMBOL = TransGeneID(gid, "ENTREZID", "SYMBOL", organism)
      paste(SYMBOL, collapse = "/")
    })
    enrichedRes@result$geneName = unlist(geneName)
  }
  return(enrichedRes)
}

