enrich.HGT = function(gene, universe, type="KEGG", organism='hsa', pvalueCutoff = 0.05,
                      pAdjustMethod = "none", minGSSize = 2, maxGSSize = 500){
  # download Kegg data
  tmp1 <- paste0("pathways_", organism," (KeGG)")
  if(!(tmp1 %in% data(package="MAGeCKFlute")$results[,"Item"])){
    gene2path_hsa=fread("http://rest.kegg.jp/link/pathway/hsa",header = F)
    names(gene2path_hsa)=c("EntrezID","PathwayID")
    pathways_hsa=fread("http://rest.kegg.jp/list/pathway/hsa",header = F)
    names(pathways_hsa)=c("PathwayID","PathwayName")
    gene2path_hsa$EntrezID=gsub("hsa:","",gene2path_hsa$EntrezID)
    gene2path_hsa$PathwayID=gsub("path:","",gene2path_hsa$PathwayID)
    pathways_hsa$PathwayID=gsub("path:","",pathways_hsa$PathwayID)
    pathways_hsa$PathwayName=gsub(" - Homo sapiens .(human.)", "", pathways_hsa$PathwayName)
    #========
    gene2path_mmu=fread("http://rest.kegg.jp/link/pathway/mmu", header = F)
    names(gene2path_mmu)=c("EntrezID","PathwayID")
    pathways_mmu=fread("http://rest.kegg.jp/list/pathway/mmu",header = F)
    names(pathways_mmu)=c("PathwayID","PathwayName")
    gene2path_mmu$EntrezID=gsub("mmu:","",gene2path_mmu$EntrezID)
    gene2path_mmu$PathwayID=gsub("path:","",gene2path_mmu$PathwayID)
    pathways_mmu$PathwayID=gsub("path:","",pathways_mmu$PathwayID)
    pathways_mmu$PathwayName=gsub(" - Mus musculus .(mouse.)", "", pathways_mmu$PathwayName)

    save(gene2path_hsa,pathways_hsa,gene2path_mmu,pathways_mmu,file="KeGG.RData")
  }

  #==========
  if(organism == "hsa"){
    gene2path=gene2path_hsa
    pathways=pathways_hsa}
  if(organism == "mmu"){
    gene2path=gene2path_mmu
    pathways=pathways_mmu}
  gene2path=as.data.frame(gene2path)
  pathways=as.data.frame(pathways)
  #============
  loginfo(paste('Running KEGG patwhay for list of entrezIDs'))
  if(type=="KEGG"){
    gene=unique(gene)
    universe=unique(universe)
    m=length(gene)
    n=length(universe)-m
    res=data.frame(ID=c(),Description=c(),GeneRatio=c(),BgRatio=c(),pvalue=c(),
                   nLogpvalue=c(),geneID=c(),geneName=c(), Count=c())
    kk=1
    for(c in pathways$PathwayID){
      pathwayEntrezID=gene2path$EntrezID[gene2path$PathwayID%in%c]
      idx1=universe %in% pathwayEntrezID
      k=length(universe[idx1])
      idx2=universe %in% gene
      q=length(universe[idx1&idx2])
      geneID=paste(universe[idx1&idx2],collapse = "/")

      if(k>=minGSSize&k<=maxGSSize&q>0){
        pvalue=phyper(q,m,n,k,lower.tail = F)
        res[kk,"ID"]=c
        res[kk,"Description"]=pathways$PathwayName[pathways$PathwayID %in% c]
        res[kk,"GeneRatio"]=paste(q,length(which(idx1)),sep="/")
        res[kk,"BgRatio"]=paste(length(which(idx1)),length(gene2path$EntrezID),sep="/")
        res[kk,"pvalue"]=pvalue
        res[kk,"nLogpvalue"]=-log10(pvalue)
        res[kk,"geneID"]=geneID

        SYMBOL = TransGeneID(universe[idx1&idx2], "ENTREZID", "SYMBOL", organism)
        geneName=paste(SYMBOL, collapse = "/")
        res[kk,"geneName"]=geneName
        res[kk,"Count"]=q

        kk=kk+1
      }
    }
    res$p.adjust=p.adjust(res$pvalue,pAdjustMethod)
    idx = which(res$pvalue<=pvalueCutoff & res$p.adjust<=pvalueCutoff)
    if(length(idx)>0){
      res = res[idx, ]
      res = res[order(res$p.adjust),]
    }else{
      res=data.frame(ID=c(),Description=c(),GeneRatio=c(),BgRatio=c(),pvalue=c(),
                     nLogpvalue=c(),geneID=c(),geneName=c(), Count=c())
    }

  }
  new("enrichResult",
      result         = res,
      pvalueCutoff   = pvalueCutoff,
      pAdjustMethod  = pAdjustMethod,
      organism       = organism,
      ontology       = type, ## as.character(x$Category[1]),
      gene           = as.character(gene),
      keytype        = "ENTREZID")
}
