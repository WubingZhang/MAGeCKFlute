#====enrichment analysis===================================
enrichment_analysis = function(geneList, universe=NULL, method=1, type="KEGG", organism="hsa",
                               pvalueCutoff = 0.05, qvalueCutoff = 1, pAdjustMethod = "BH",
                               minGSSize = 2, maxGSSize = 500, plotTitle=NULL,gridColour="blue"){

  result=list()
  type=toupper(type[1])
  methods = c("ORT", "GSEA", "DAVID", "GOstats", "HyperGeometric")
  names(methods) = toupper(methods)
  if(class(method)=="character"){method = toupper(method)}
  method = methods[method]
  #======================================================================================
  p1=ggplot()
  p1=p1+geom_text(aes(x=0,y=0,label="Less than 10 genes"),size=6)
  p1=p1+labs(title=plotTitle)
  p1=p1+theme(plot.title = element_text(size=10))
  p1=p1+theme_void()
  p1=p1+theme(plot.title = element_text(hjust = 0.5))
  if(length(geneList)<10){
    result$enrichRes = NULL
    result$gridPlot = p1
    return(result)
  }
  #====Gene Set Enrichment Analysis=======================================================
  if(method == "GSEA"){
    geneList=geneList[order(geneList,decreasing = T)]
    enrichRes <- enrich.GSE(geneList, type = type, pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod,
                            organism=organism, minGSSize = minGSSize, maxGSSize = maxGSSize)
    result$enrichRes = enrichRes
    gridPlot <- enrichment_plot_GSE(enrichRes@result, plotTitle, gridColour=gridColour)
    result$gridPlot = gridPlot
    return(result)
  }
  #====Over-Representation Analysis======================================================
  if(method == "ORT"){
    enrichRes <- enrich.ORT(gene = geneList, universe = universe, type = type, organism=organism,
                            pvalueCutoff=pvalueCutoff, qvalueCutoff = qvalueCutoff, pAdjustMethod = pAdjustMethod,
                            minGSSize = minGSSize, maxGSSize = maxGSSize)
  }
  #=============DAVID=======================================================================
  if(method == "DAVID"){
    if(type == "KEGG"){
      annotation = "KEGG_PATHWAY"
    }else if(type %in% c("BP", "CC", "MF")){
      annotation = paste("GOTERM", type, "FAT", sep="_")
    }else if(type == "DO"){
      annotation = "OMIM_DISEASE"
    }else{annotation = type}

    enrichRes = enrich.DAVID(gene = geneList, universe = universe,
                             minGSSize = minGSSize, maxGSSize = maxGSSize, annotation  = annotation,
                             pvalueCutoff  = pvalueCutoff, pAdjustMethod = pAdjustMethod, qvalueCutoff= qvalueCutoff)
  }
  #===============================GOstats enrichment analysis============================
  if(method == "GOstats"){
    enrichRes = enrich.GOstats(gene = geneList, universe = universe, type  = type, organism=organism,
                             pvalueCutoff  = pvalueCutoff, pAdjustMethod = pAdjustMethod)
  }
  #==================================HyperGeometric test=================================
  if(method == "HyperGeometric"){
    enrichRes = enrich.HGT(gene = geneList, universe = universe, type  = type, organism=organism,
                             pvalueCutoff  = pvalueCutoff, pAdjustMethod = pAdjustMethod,
                           minGSSize = minGSSize, maxGSSize = maxGSSize)
  }
  result$enrichRes = enrichRes
  gridPlot <- enrichment_plot(enrichRes@result, plotTitle, gridColour=gridColour)
  result$gridPlot = gridPlot

  return(result)
}



