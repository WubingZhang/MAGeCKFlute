#====enrichment analysis===================================
enrichment_analysis = function(geneList, genes, method=1, type=c("KEGG", "BP", "MF", "CC", "DO", "NCG", "MKEGG"),
                               pvalueCutoff = 0.05, qvalueCutoff = 1, pAdjustMethod = "BH",
                               minGSSize = 2, maxGSSize = 500, plotTitle=NULL,gridColour="blue"){

  result=list()
  type=toupper(type[1])
  geneList=geneList[order(geneList,decreasing = T)]
  methods = c("ORT", "GSEA", "DAVID", "GOstats", "HyperGeometric")
  names(methods) = toupper(methods)
  if(class(method)=="character"){method = toupper(method)}
  method = methods[method]
  #====Over-Representation Analysis===============
  if(method == "ORT"){
    enrichRes <- enrich.ORT(gene = genes, universe = names(geneList), type = type,
                            pvalueCutoff=pvalueCutoff, qvalueCutoff = qvalueCutoff, pAdjustMethod = pAdjustMethod,
                            minGSSize = minGSSize, maxGSSize = maxGSSize)
    result$enrichRes = enrichRes
    gridPlot <- enrichment_plot(enrichRes@result, plotTitle, gridColour=gridColour)
    result$gridPlot = gridPlot
  }
  #====Gene Set Enrichment Analysis===============
  if(method == "GSEA"){
    enrichRes <- enrich.GSE(geneList, type = type, pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod,
                            minGSSize = minGSSize, maxGSSize = maxGSSize)
    result$enrichRes = enrichRes
    gridPlot <- enrichment_plot_GSE(enrichRes@result, plotTitle, gridColour=gridColour)
    result$gridPlot = gridPlot
  }
  #====DAVID===============
  if(method == "DAVID"){
    if(type == "KEGG"){
      annotation = "KEGG_PATHWAY"
    }else if(type %in% c("BP", "CC", "MF")){
      annotation = paste("GOTERM", type, "FAT", sep="_")
    }else if(type == "DO"){
      annotation = "OMIM_DISEASE"
    }else{annotation = type}

    enrichRes = enrich.DAVID(gene = genes, universe = names(geneList),
                             minGSSize = minGSSize, maxGSSize = maxGSSize, annotation  = annotation,
                             pvalueCutoff  = pvalueCutoff, pAdjustMethod = pAdjustMethod, qvalueCutoff= qvalueCutoff)
    result$enrichRes = enrichRes
    gridPlot <- enrichment_plot(enrichRes@result, plotTitle, gridColour=gridColour)
    result$gridPlot = gridPlot
  }
  #===GOstats enrichment analysis==============
  if(method == "GOstats"){
    enrichRes = enrich.GOstats(gene = genes, universe = names(geneList), type  = type,
                             pvalueCutoff  = pvalueCutoff, pAdjustMethod = pAdjustMethod)
    result$enrichRes = enrichRes
    gridPlot <- enrichment_plot(enrichRes@result, plotTitle, gridColour=gridColour)
    result$gridPlot = gridPlot
  }
  #======HyperGeometric test==============
  if(method == "HyperGeometric"){
    enrichRes = enrich.HGT(gene = genes, universe = names(geneList), type  = type,
                             pvalueCutoff  = pvalueCutoff, pAdjustMethod = pAdjustMethod)
    result$enrichRes = enrichRes
    gridPlot <- enrichment_plot(enrichRes@result, plotTitle, gridColour=gridColour)
    result$gridPlot = gridPlot
  }

  return(result)
}

enrich.ORT <- function(gene, universe=NULL, type=c("KEGG", "BP", "MF", "CC", "DO", "NCG", "MKEGG"),
                       pvalueCutoff = 0.05, qvalueCutoff = 0.2, pAdjustMethod = "BH",
                       minGSSize = 2, maxGSSize = 500){
  loginfo(paste('Running Over-Representation Test for list of entrezIDs'))
  if(type %in% c("BP", "CC", "MF")){
    enrichedRes = enrichGO(gene = gene, universe = universe, ont = type, OrgDb=org.Hs.eg.db,
              pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff,
              minGSSize=minGSSize, maxGSSize=maxGSSize)
    enrichedRes = simplify(enrichedRes, cutoff=0.7, by="p.adjust", select_fun=min)
  }
  if(type == "KEGG"){
    enrichedRes = enrichKEGG(gene=gene,  universe=universe,
                           pAdjustMethod = pAdjustMethod,pvalueCutoff = pvalueCutoff,qvalueCutoff=qvalueCutoff,
                           minGSSize = minGSSize, maxGSSize = maxGSSize, use_internal_data=F)
  }
  if(type == "DO"){
    enrichedRes = enrichDO(gene=gene,  universe=universe,
                           pAdjustMethod = pAdjustMethod,pvalueCutoff = pvalueCutoff,qvalueCutoff=qvalueCutoff,
                           minGSSize = minGSSize, maxGSSize = maxGSSize)
  }
  if(type == "MKEGG"){
    enrichedRes = enrichMKEGG(gene=gene,  universe=universe,
                           pAdjustMethod = pAdjustMethod,pvalueCutoff = pvalueCutoff,qvalueCutoff=qvalueCutoff,
                           minGSSize = minGSSize, maxGSSize = maxGSSize)
  }
  if(type == "NCG"){  loginfo(paste('Running KEGG for list of entrezIDs'))
    enrichedRes = enrichNCG(gene=gene,  universe=universe,
                              pAdjustMethod = pAdjustMethod,pvalueCutoff = pvalueCutoff,qvalueCutoff=qvalueCutoff,
                              minGSSize = minGSSize, maxGSSize = maxGSSize)
  }
  return(enrichedRes)
}

enrich.GSE <- function(geneList, type=c("KEGG", "BP", "MF", "CC", "DO", "NCG", "MKEGG"),
                       minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "BH"){
  loginfo(paste('Running GSEA for list of entrezIDs'))
  #geneList:	order ranked geneList
  if(type %in% c("BP", "CC", "MF")){
    enrichedRes = gseGO(geneList=geneList, ont = DS, OrgDb=org.Hs.eg.db,
             minGSSize = minGSSize, maxGSSize = maxGSSize,
             pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod)
  }
  if(type == "KEGG"){
    # download Kegg data
    # KEGGRdata="KeGG.RData"

    # if(file.exists(KEGGRdata)){
    #   load(KEGGRdata)
    # }else{
    #   gene2path=fread(paste("http://rest.kegg.jp/link/pathway/",organism,sep=""),header = F)
    #   names(gene2path)=c("EntrezID","PathwayID")
    #   pathways=fread(paste("http://rest.kegg.jp/list/pathway/",organism,sep=""),header = F)
    #   names(pathways)=c("PathwayID","PathwayName")
    #   save(gene2path,pathways,file =KEGGRdata)
    # }
    data("KeGG")

    gene2path$EntrezID=gsub("hsa:","",gene2path$EntrezID)
    gene2path$PathwayID=gsub("path:","",gene2path$PathwayID)
    pathways$PathwayID=gsub("path:","",pathways$PathwayID)
    pathways$PathwayName=gsub(" - Homo sapiens .(human.)", "", pathways$PathwayName)
    pathways$PathwayName=gsub(" - Mus musculus .(mouse.)", "", pathways$PathwayName)
    gene2path=as.data.frame(gene2path)
    pathways=as.data.frame(pathways)

    enrichedRes = GSEA(geneList=geneList, minGSSize = minGSSize, maxGSSize = maxGSSize,
                        pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod,
                       TERM2GENE=gene2path[,c("PathwayID","EntrezID")], TERM2NAME=pathways)
  }
  if(type == "DO"){
    enrichedRes = gseDO(geneList=geneList, minGSSize = minGSSize, maxGSSize = maxGSSize,
                          pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod)
  }
  if(type == "MKEGG"){
    enrichedRes = gseMKEGG(geneList=geneList, minGSSize = minGSSize, maxGSSize = maxGSSize,
                          pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod)
  }
  if(type == "NCG"){
    enrichedRes = gseMKEGG(geneList=geneList, minGSSize = minGSSize, maxGSSize = maxGSSize,
                           pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod)
  }
  return(enrichedRes)
}


enrich.DAVID <- function(gene, universe=NULL, david.user="ma.tongji@gmail.com",
                        minGSSize = 2, maxGSSize = 500, annotation  = "GOTERM_BP_FAT",
                        pvalueCutoff  = 0.05, pAdjustMethod = "BH", qvalueCutoff= 0.2){

  # A job with more than 3000 genes to generate gene or term cluster report will not be handled by DAVID due to resource limit.
  # No more than 200 jobs in a day from one user or computer.
  # DAVID Team reserves right to suspend any improper uses of the web service without notice.
  ##' enrichment analysis by DAVID

  loginfo(paste('Running DAVID for list of entrezIDs'))

  Count <- List.Total <- Pop.Hits <- Pop.Total <- NULL

  pAdjustMethod <- match.arg(pAdjustMethod, c("bonferroni", "BH"))

  david.pkg <- "RDAVIDWebService"
  pkgs <- installed.packages()[,1]
  if (! david.pkg %in% pkgs) {
    stop("You should have RDAVIDWebService package installed before using enrichDAVID...")
  }

  require(david.pkg, character.only=TRUE)
  DAVIDWebService <- eval(parse(text="DAVIDWebService"))
  addList <- eval(parse(text="addList"))
  setAnnotationCategories <- eval(parse(text="setAnnotationCategories"))
  getFunctionalAnnotationChart <- eval(parse(text="getFunctionalAnnotationChart"))
  getSpecieNames <- eval(parse(text="getSpecieNames"))
  getIdTypes <- eval(parse(text="getIdTypes"))

  david <- DAVIDWebService$new(email=david.user,
                               url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

  ## addList will throw error if idType is not match.
  ## use match.arg to check before addList make it more readable

  idType <- match.arg("ENTREZ_GENE_ID", getIdTypes(david))

  david.res <- addList(david, gene, idType=idType, listName="clusterProfiler", listType="Gene")
  if(!is.null(universe)){
    david.res <- addList(david, universe, idType=idType, listName="User_assign", listType="Background")
  }

  if (david.res$inDavid == 0) {
    stop("All id can not be mapped. Please check 'idType' parameter...")
  }

  setAnnotationCategories(david, annotation)
  x <- getFunctionalAnnotationChart(david, threshold=1, count=minGSSize)

  if (length(x@.Data) == 0) {
    warning("No significant enrichment found...")
    return(NULL)
  }

  term <- x$Term
  if (length(grep("~", term[1])) == 0) {
    sep <- ":"
  } else {
    sep <- "~"
  }
  term.list <- sapply(term, function(y) strsplit(y, split=sep))
  term.df <- do.call("rbind", term.list)
  ID <- term.df[,1]
  Description <- term.df[,2]
  GeneRatio <- with(x, paste(Count, List.Total, sep="/"))
  BgRatio <- with(x, paste(Pop.Hits, Pop.Total, sep="/"))
  Over <- data.frame(ID          = ID,
                     Description = Description,
                     GeneRatio   = GeneRatio,
                     BgRatio     = BgRatio,
                     pvalue      = x$PValue,
                     stringsAsFactors = FALSE)
  row.names(Over) <- ID

  if (pAdjustMethod == "bonferroni") {
    Over$p.adjust <- x$Bonferroni
  } else {
    Over$p.adjust <- x$Benjamini
  }

  qobj <- tryCatch(qvalue(p=Over$pvalue, lambda=0.05, pi0.method="bootstrap"),
                   error=function(e) NULL)
  if (class(qobj) == "qvalue") {
    qvalues <- qobj$qvalues
  } else {
    qvalues <- NA
  }
  Over$qvalue <- qvalues
  Over$geneID <- gsub(",\\s*", "/", x$Genes)
  Over$Count <- x$Count

  Over <- Over[ Over$pvalue <= pvalueCutoff, ]
  Over <- Over[ Over$p.adjust <= pvalueCutoff, ]
  if (! any(is.na(Over$qvalue))) {
    Over <- Over[ Over$qvalue <= qvalueCutoff, ]
  }

  org <- getSpecieNames(david)
  org <- gsub("\\(.*\\)", "", org)

  ## gc <- strsplit(Over$geneID, "/")
  ## names(gc) <- Over$ID

  if (!is.na(maxGSSize) || !is.null(maxGSSize)) {
    idx <- as.numeric(sub("/\\d+", "", Over$BgRatio)) <= maxGSSize
    Over <- Over[idx,]
  }

  new("enrichResult",
      result         = Over,
      pvalueCutoff   = pvalueCutoff,
      pAdjustMethod  = pAdjustMethod,
      organism       = org,
      ontology       = annotation, ## as.character(x$Category[1]),
      gene           = as.character(gene),
      keytype        = idType)
}

enrich.GOstats <- function(gene, universe=NULL, type=c("KEGG", "BP", "MF", "CC"), pvalueCutoff = 0.05, pAdjustMethod = "BH"){
  #======
  DS = toupper(type[1])
  over.sum=NULL
  if(DS == "KEGG"){
    loginfo(paste('Running KEGG for list of entrezIDs'))
    params <- new("KEGGHyperGParams", categoryName="KEGG",
                  geneIds=gene, universeGeneIds=universe,
                  annotation="org.Hs.eg.db", pvalueCutoff=pvalueCutoff,testDirection="over")
    #loginfo('    Starting HyperG test')
    over <- hyperGTest(params)
    #loginfo('    HyperG test done')
    over.sum <- summary(over)

    glist <- geneIdsByCategory(over)
    glist <- sapply(glist, function(.ids) {
      .sym <- .ids
      # .sym <- mget(.ids, envir=org.Hs.egSYMBOL, ifnotfound=NA)
      .sym[is.na(.sym)] <- .ids[is.na(.sym)]
      paste(.sym, collapse="/")
    })
    over.sum$geneID <- glist[as.character(over.sum$KEGGID)]
    # over.sum$Symbols <- glist[as.character(over.sum$KEGGID)]
  }
  if(DS %in% c("BP", "CC", "MF")){
    # gene ontology background
    loginfo(paste('Running GO', DS, 'for list of entrezIDs'))
    params <- new("GOHyperGParams", annotation="org.Hs.eg.db",
                  geneIds = gene, universeGeneIds = universe,
                  ontology=DS, pvalueCutoff=pvalueCutoff, testDirection="over")
    over <- hyperGTest(params)
    over.sum <- summary(over)
    glist <- geneIdsByCategory(over)
    glist <- sapply(glist, function(.ids) {
      .sym <- .ids
      # .sym <- mget(.ids, envir=org.Hs.egSYMBOL, ifnotfound=NA)
      .sym[is.na(.sym)] <- .ids[is.na(.sym)]
      paste(.sym, collapse="/")
    })
    over.sum$geneID <- glist[as.character(over.sum[,1])]
  }
  result = list()
  result$ID = over.sum[,1]
  result$Description = over.sum$Term
  result$GeneRatio = paste0(over.sum$Count, "/", over.sum$Size)
  result$BgRatio = paste0(over.sum$Size, "/", length(universe))
  result$pvalue = over.sum$Pvalue
  result$p.adjust = p.adjust(result$pvalue, method=pAdjustMethod)
  result$geneID = over.sum$geneID
  result$Count = over.sum$Count
  result = as.data.frame(result)
  idx = which(result$pvalue<=pvalueCutoff & result$p.adjust<=pvalueCutoff)
  result = result[idx, ]
  result = result[order(result$p.adjust),]
  new("enrichResult",
      result         = result,
      pvalueCutoff   = pvalueCutoff,
      pAdjustMethod  = pAdjustMethod,
      organism       = "Homo sapiens",
      ontology       = type, ## as.character(x$Category[1]),
      gene           = as.character(gene),
      keytype        = "ENTREZID")
}

enrich.HGT = function(gene, universe, type="KEGG", organism='hsa', pvalueCutoff = 0.05, pAdjustMethod = "BH"){
  # download Kegg data
  # KEGGRdata=file.path(PATH, "../data/KeGG.RData")
  # cat(KEGGRdata)
  # if(file.exists(KEGGRdata)){
  #   load(KEGGRdata)
  # }else{
  #   gene2path=fread(paste("http://rest.kegg.jp/link/pathway/",organism,sep=""),header = F)
  #   names(gene2path)=c("EntrezID","PathwayID")
  #   pathways=fread(paste("http://rest.kegg.jp/list/pathway/",organism,sep=""),header = F)
  #   names(pathways)=c("PathwayID","PathwayName")
  #   save(gene2path,pathways,file =KEGGRdata)
  # }
  data("KeGG")
  gene2path$EntrezID=gsub("hsa:","",gene2path$EntrezID)
  gene2path$PathwayID=gsub("path:","",gene2path$PathwayID)
  pathways$PathwayID=gsub("path:","",pathways$PathwayID)
  pathways$PathwayName=gsub(" - Homo sapiens .(human.)", "", pathways$PathwayName)
  pathways$PathwayName=gsub(" - Mus musculus .(mouse.)", "", pathways$PathwayName)
  gene2path=as.data.frame(gene2path)
  pathways=as.data.frame(pathways)
  #============
  loginfo(paste('Running KEGG patwhay for list of entrezIDs'))
  res=NULL
  if(type=="KEGG"){
    m=length(gene)
    n=length(universe)-m
    res=data.frame(ID=c(),Description=c(),GeneRatio=c(),BgRatio=c(),pvalue=c(),nLogpvalue=c(),geneID=c(),Count=c())
    kk=1
    for(c in pathways$PathwayID){
      pathwayEntrezID=gene2path$EntrezID[gene2path$PathwayID%in%c]
      idx1=universe %in% pathwayEntrezID
      k=length(universe[idx1])
      idx2=universe %in% gene
      q=length(universe[idx1&idx2])
      geneID=paste(universe[idx1&idx2],collapse = "/")

      if(k>10&q>0){
        pvalue=phyper(q,m,n,k,lower.tail = F)
        res[kk,"ID"]=c
        res[kk,"Description"]=pathways$PathwayName[pathways$PathwayID %in% c]
        res[kk,"GeneRatio"]=paste(q,length(which(idx1)),sep="/")
        res[kk,"BgRatio"]=paste(length(which(idx1)),length(gene2path$EntrezID),sep="/")
        res[kk,"pvalue"]=pvalue
        res[kk,"nLogpvalue"]=-log10(pvalue)
        res[kk,"geneID"]=geneID
        res[kk,"Count"]=q

        kk=kk+1

      }
    }
    res$p.adjust=p.adjust(res$pvalue,pAdjustMethod)
    idx = which(res$pvalue<=pvalueCutoff & res$p.adjust<=pvalueCutoff)
    res = res[idx, ]
    res = res[order(res$p.adjust),]
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

enrichment_plot=function(enrichment,plotTitle=NULL,gridColour="blue",termNum=10,charLength=50){
  #The column of Description, ID, p.adjust, and Count are neccessary.
  enrichment$logP = -log10(enrichment$p.adjust)
  enrichment = enrichment[!is.na(enrichment$ID),]
  enrichment=enrichment[!duplicated(enrichment$Description),]
  enrichment = enrichment[order(enrichment$logP,decreasing = T), ]

  if(nrow(enrichment)==0){
    p1=ggplot()
    p1=p1+geom_text(aes(x=0,y=0,label="No enriched terms"),size=6)
    p1=p1+labs(title=plotTitle)
    p1=p1+theme(plot.title = element_text(size=10))
    p1=p1+theme_void()
    p1=p1+theme(plot.title = element_text(hjust = 0.5))
    return(p1)
  }

  if(nrow(enrichment)>=termNum){
    enrichment=enrichment[1:termNum,]
  }
  #normalize term description
  {
    terms=as.character(enrichment$Description)
    terms=lapply(terms,function(x,k){
      x=as.character(x)
      x=paste(toupper(substr(x,1,1)),substr(x,2,nchar(x)),sep="")
      if(nchar(x)>k){
        x=substr(x,start=1,stop=k)
        x=gsub("(\\w+)$", "...", x)
      }
      return(x)
    },charLength)
    enrichment$Description=do.call(rbind,terms)
  }
  idx = duplicated(enrichment$Description)
  enrichment=enrichment[!idx,]
  enrichment$Name=factor(enrichment$Description,levels=rev(enrichment$Description))

  p1=ggplot(data=enrichment,aes(x=logP,y=Name,size=Count,colour = factor(Count)))
  p1=p1+geom_point()
  p1=p1+guides(color=F)
  p1 <- p1+theme(panel.grid.major=element_line(colour=gridColour),
                 panel.grid.minor=element_blank(),
                 panel.background=element_blank())
  p1 <- p1+xlab("-log10(Adjust.pvalue)")+ylab("")
  p1 <- p1+labs(title=plotTitle)
  p1 <- p1+theme(axis.text.x=element_text(size=6, face="plain", colour='black'))
  p1 <- p1+theme(axis.text.y=element_text(size=6, face="plain", colour='black'))
  p1=p1+theme(legend.position="bottom")+theme(plot.title = element_text(hjust = 0.5,size=10))
  #p1
  return(p1)
}

enrichment_plot_GSE=function(enrichment,plotTitle=NULL,gridColour="blue",termNum=10,charLength=50){
  #The column of Description, ID, p.adjust, and Count are neccessary.
  enrichment$logP = -log10(enrichment$p.adjust)
  enrichment = enrichment[!is.na(enrichment$ID),]
  enrichment=enrichment[!duplicated(enrichment$Description),]
  enrichment = enrichment[order(enrichment$logP,decreasing = T), ]

  if(nrow(enrichment)==0){
    p1=ggplot()
    p1=p1+geom_text(aes(x=0,y=0,label="No enriched terms"),size=6)
    p1=p1+labs(title=plotTitle)
    p1=p1+theme(plot.title = element_text(size=10))
    p1=p1+theme_void()
    p1=p1+theme(plot.title = element_text(hjust = 0.5))
    return(p1)
  }

  if(nrow(enrichment)>=termNum){
    enrichment=enrichment[1:termNum,]
  }
  #normalize term description
  {
    terms=as.character(enrichment$Description)
    terms=lapply(terms,function(x,k){
      x=as.character(x)
      x=paste(toupper(substr(x,1,1)),substr(x,2,nchar(x)),sep="")
      if(nchar(x)>k){
        x=substr(x,start=1,stop=k)
        x=gsub("(\\w+)$", "...", x)
      }
      return(x)
    },charLength)
    enrichment$Description=do.call(rbind,terms)
  }
  idx = duplicated(enrichment$Description)
  enrichment=enrichment[!idx,]
  enrichment$Name=factor(enrichment$Description,levels=rev(enrichment$Description))

  p1=ggplot(data=enrichment,aes(x=logP,y=Name,size=NES,colour = factor(NES)))
  p1=p1+geom_point()
  p1=p1+guides(color=F)
  p1 <- p1+theme(panel.grid.major=element_line(colour=gridColour),
                 panel.grid.minor=element_blank(),
                 panel.background=element_blank())
  p1 <- p1+xlab("-log10(Adjust.pvalue)")+ylab("")
  p1 <- p1+labs(title=plotTitle)
  p1 <- p1+theme(axis.text.x=element_text(size=6, face="plain", colour='black'))
  p1 <- p1+theme(axis.text.y=element_text(size=6, face="plain", colour='black'))
  p1=p1+theme(legend.position="bottom")+theme(plot.title = element_text(hjust = 0.5,size=10))
  #p1
  return(p1)
}
