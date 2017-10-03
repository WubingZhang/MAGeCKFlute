
enrich.GSE <- function(geneList, type= "KEGG", organism='hsa', minGSSize = 10, maxGSSize = 500,
                       pvalueCutoff = 0.05, pAdjustMethod = "none"){

  geneList = sort(geneList, decreasing = T)
  loginfo(paste('Running GSEA for list of entrezIDs'))
  #geneList:	order ranked geneList
  if(type == "KEGG"){
    # download Kegg data
    tmp1 <- paste0("pathways_", organism," (KeGG)")
    if(!(tmp1 %in% data(package="MAGeCKFlute")$results[,"Item"])){
      gene2path_hsa=fread("http://rest.kegg.jp/link/pathway/hsa",header = F)
      names(gene2path_hsa)=c("EntrezID","PathwayID")
      pathways_hsa=fread("http://rest.kegg.jp/list/pathway/hsa",header = F)
      names(pathways_hsa)=c("PathwayID","PathwayName")
      gene2path_hsa$PathwayID=gsub("path:","",gene2path_hsa$PathwayID)
      pathways_hsa$PathwayID=gsub("path:","",pathways_hsa$PathwayID)
      pathways_hsa$PathwayName=gsub(" - Homo sapiens .(human.)", "", pathways_hsa$PathwayName)

      #========
      gene2path_mmu=fread("http://rest.kegg.jp/link/pathway/mmu",header = F)
      names(gene2path_mmu)=c("EntrezID","PathwayID")
      pathways_mmu=fread("http://rest.kegg.jp/list/pathway/mmu",header = F)
      names(pathways_mmu)=c("PathwayID","PathwayName")
      gene2path_mmu$PathwayID=gsub("path:","",gene2path_mmu$PathwayID)
      pathways_mmu$PathwayID=gsub("path:","",pathways_mmu$PathwayID)
      pathways_mmu$PathwayName=gsub(" - Mus musculus .(mouse.)", "", pathways_mmu$PathwayName)

      save(gene2path_hsa,pathways_hsa,gene2path_mmu,pathways_mmu,file="KeGG.RData")
    }

    #==================
    if(organism == "hsa"){
      gene2path=gene2path_hsa
      pathways=pathways_hsa}
    if(organism == "mmu"){
      gene2path=gene2path_mmu
      pathways=pathways_mmu}

    enrichedRes = GSEA(geneList=geneList, minGSSize = minGSSize, maxGSSize = maxGSSize,
                       pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod,
                       TERM2GENE=gene2path[,c("PathwayID","EntrezID")], TERM2NAME=pathways)
  }

  if(type %in% c("BP", "CC", "MF")){
    orgdb = c("org.Hs.eg.db", "org.Mm.eg.db")
    names(orgdb) = c("hsa", "mmu")
    enrichedRes = gseGO(geneList=geneList, ont = type, OrgDb=orgdb[organism],
                        minGSSize = minGSSize, maxGSSize = maxGSSize,
                        pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod)
  }

  if(type == "DO"){
    enrichedRes = gseDO(geneList=geneList, minGSSize = minGSSize, maxGSSize = maxGSSize,
                        pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod)
  }
  if(type == "MKEGG"){
    enrichedRes = gseMKEGG(geneList=geneList, organism = organism, minGSSize = minGSSize, maxGSSize = maxGSSize,
                           pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod)
  }
  if(type == "NCG"){
    enrichedRes = gseNCG(geneList=geneList, minGSSize = minGSSize, maxGSSize = maxGSSize,
                           pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod)
  }
  if(!is.null(enrichedRes)){
    geneID = strsplit(enrichedRes@result$core_enrichment, "/")
    geneName = lapply(geneID, function(gid){
      SYMBOL = TransGeneID(gid, "ENTREZID", "SYMBOL", organism)
      paste(SYMBOL, collapse = "/")
    })
    enrichedRes@result$geneName = unlist(geneName)
  }

  return(enrichedRes)
}

