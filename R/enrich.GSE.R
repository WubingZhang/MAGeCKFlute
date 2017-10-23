
enrich.GSE <- function(geneList, type= "KEGG", organism='hsa', minGSSize = 10, maxGSSize = 500,
                       pvalueCutoff = 0.05, pAdjustMethod = "BH"){

  geneList = sort(geneList, decreasing = TRUE)
  loginfo(paste('Running GSEA for list of entrezIDs'))
  #geneList:	order ranked geneList
  if(type == "KEGG"){
    # download Kegg data
    organism = getOrg(organism, onlyLib = TRUE)$organism
    tmp1 <- paste0("pathways_", organism," (KeGG)")
    if(!(tmp1 %in% data(package="MAGeCKFlute")$results[,"Item"])){
      gene2path=fread(paste0("http://rest.kegg.jp/link/pathway/",organism),header = FALSE)
      names(gene2path)=c("EntrezID","PathwayID")
      pathways=fread(paste0("http://rest.kegg.jp/list/pathway/",organism),header = FALSE)
      names(pathways)=c("PathwayID","PathwayName")
      gene2path$PathwayID=gsub("path:","",gene2path$PathwayID)
      pathways$PathwayID=gsub("path:","",pathways$PathwayID)
      pathways$PathwayName=gsub(" - Homo sapiens .(human.)", "", pathways$PathwayName)

      # #========
      # gene2path_mmu=fread("http://rest.kegg.jp/link/pathway/mmu",header = FALSE)
      # names(gene2path_mmu)=c("EntrezID","PathwayID")
      # pathways_mmu=fread("http://rest.kegg.jp/list/pathway/mmu",header = FALSE)
      # names(pathways_mmu)=c("PathwayID","PathwayName")
      # gene2path_mmu$PathwayID=gsub("path:","",gene2path_mmu$PathwayID)
      # pathways_mmu$PathwayID=gsub("path:","",pathways_mmu$PathwayID)
      # pathways_mmu$PathwayName=gsub(" - Mus musculus .(mouse.)", "", pathways_mmu$PathwayName)
      #
      # save(gene2path_hsa,pathways_hsa,gene2path_mmu,pathways_mmu,file="KeGG.RData")
    }else if(organism == "hsa"){
      gene2path=gene2path_hsa
      pathways=pathways_hsa
    }else if(organism == "mmu"){
      gene2path=gene2path_mmu
      pathways=pathways_mmu}

    enrichedRes = GSEA(geneList=geneList, minGSSize = minGSSize, maxGSSize = maxGSSize,
                       pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod,
                       TERM2GENE=gene2path[,c("PathwayID","EntrezID")], TERM2NAME=pathways)
  }

  if(type %in% c("BP", "CC", "MF")){
    orgdb = getOrg(organism, onlyLib = TRUE)$pkg
    enrichedRes = gseGO(geneList=geneList, ont = type, OrgDb=orgdb,
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
  if(!is.null(enrichedRes) && nrow(enrichedRes@result)>0){
    geneID = strsplit(enrichedRes@result$core_enrichment, "/")
    geneName = lapply(geneID, function(gid){
      SYMBOL = TransGeneID(gid, "ENTREZID", "SYMBOL", organism)
      paste(SYMBOL, collapse = "/")
    })
    enrichedRes@result$geneName = unlist(geneName)
  }

  return(enrichedRes)
}

