#enrichment for square grouped genes
EnrichSquare <- function(beta, ctrlname="Control",treatname="Treatment", pvalue=0.05,enrich_method="Hypergeometric",
                         organism="hsa", adjust="none", filename=NULL, out.dir="."){
  loginfo("Enrichment analysis of 9 Square grouped genes ...")
  gg=beta

  idx=gg$group=="Group1"
  genes=as.character(gg$ENTREZID[idx])
  universe=as.character(gg$ENTREZID)

  #====GO_KEGG_enrichment=====
  kegg1=enrichment_analysis(geneList = genes, universe=universe, method = enrich_method,type = "KEGG",organism=organism,
                            pvalueCutoff = pvalue, plotTitle="KEGG: Group1", pAdjustMethod = adjust)
  bp1=enrichment_analysis(geneList = genes, universe=universe,  method = "ORT", type = "BP",organism=organism,
                          pvalueCutoff = pvalue, plotTitle="BP: Group1", pAdjustMethod = adjust)


  idx=gg$group=="Group2"
  genes=as.character(gg$ENTREZID[idx])
  #====GO_KEGG_enrichment=====
  kegg2=enrichment_analysis(geneList = genes, universe=universe, method = enrich_method,type = "KEGG",organism=organism,
                            pvalueCutoff = pvalue, plotTitle="KEGG: Group2", pAdjustMethod = adjust)
  bp2=enrichment_analysis(geneList = genes, universe=universe, method = "ORT", type = "BP",organism=organism,
                          pvalueCutoff = pvalue, plotTitle="BP: Group2", pAdjustMethod = adjust)

  idx=gg$group=="Group3"
  genes=as.character(gg$ENTREZID[idx])
  #====GO_KEGG_enrichment=====
  kegg3=enrichment_analysis(geneList = genes, universe=universe, method = enrich_method,type = "KEGG",organism=organism,
                            pvalueCutoff = pvalue, plotTitle="KEGG: Group3", pAdjustMethod = adjust)
  bp3=enrichment_analysis(geneList = genes, universe=universe, method = "ORT", type = "BP",organism=organism,
                          pvalueCutoff = pvalue, plotTitle="BP: Group3", pAdjustMethod = adjust)

  idx=gg$group=="Group4"
  genes=as.character(gg$ENTREZID[idx])
  #====GO_KEGG_enrichment=====
  kegg4=enrichment_analysis(geneList = genes, universe=universe, method = enrich_method,type = "KEGG",organism=organism,
                            pvalueCutoff = pvalue, plotTitle="KEGG: Group4", pAdjustMethod = adjust)
  bp4=enrichment_analysis(geneList = genes, universe=universe, method = "ORT", type = "BP",organism=organism,
                          pvalueCutoff = pvalue, plotTitle="BP: Group4", pAdjustMethod = adjust)

  idx1=gg$group=="Group1"
  idx2=gg$group=="Group3"
  idx = idx1|idx2
  genes=as.character(gg$ENTREZID[idx])
  #====GO_KEGG_enrichment=====
  kegg13=enrichment_analysis(geneList = genes, universe=universe, method = enrich_method,type = "KEGG",organism=organism,
                             pvalueCutoff = pvalue, plotTitle="KEGG: Group1&3", pAdjustMethod = adjust)
  bp13=enrichment_analysis(geneList = genes, universe=universe, method = "ORT", type = "BP",organism=organism,
                           pvalueCutoff = pvalue, plotTitle="BP: Group1&3", pAdjustMethod = adjust)


  idx1=gg$group=="Group1"
  idx2=gg$group=="Group4"
  idx = idx1|idx2
  genes=as.character(gg$ENTREZID[idx])
  #====GO_KEGG_enrichment=====
  kegg14=enrichment_analysis(geneList = genes, universe=universe, method = enrich_method,type = "KEGG",organism=organism,
                             pvalueCutoff = pvalue, plotTitle="KEGG: Group1&4", pAdjustMethod = adjust)
  bp14=enrichment_analysis(geneList = genes, universe=universe, method = "ORT", type = "BP",organism=organism,
                           pvalueCutoff = pvalue, plotTitle="BP: Group1&4", pAdjustMethod = adjust)

  idx1=gg$group=="Group2"
  idx2=gg$group=="Group3"
  idx = idx1|idx2
  genes=as.character(gg$ENTREZID[idx])
  #====GO_KEGG_enrichment=====
  kegg23=enrichment_analysis(geneList = genes, universe=universe, method = enrich_method,type = "KEGG",organism=organism,
                             pvalueCutoff = pvalue, plotTitle="KEGG: Group2&3", pAdjustMethod = adjust)
  bp23=enrichment_analysis(geneList = genes, universe=universe, method = "ORT", type = "BP",organism=organism,
                           pvalueCutoff = pvalue, plotTitle="BP: Group2&3", pAdjustMethod = adjust)

  idx1=gg$group=="Group2"
  idx2=gg$group=="Group4"
  idx = idx1|idx2
  genes=as.character(gg$ENTREZID[idx])

  #====GO_KEGG_enrichment=====
  kegg24=enrichment_analysis(geneList = genes, universe=universe, method = enrich_method,type = "KEGG",organism=organism,
                             pvalueCutoff = pvalue, plotTitle="KEGG: Group2&4", pAdjustMethod = adjust)
  bp24=enrichment_analysis(geneList = genes, universe=universe, method = "ORT", type = "BP",organism=organism,
                           pvalueCutoff = pvalue, plotTitle="BP: Group2&4", pAdjustMethod = adjust)

  idx1=gg$group=="Group1"
  idx2=gg$group=="Group2"
  idx3=gg$group=="Group3"
  idx4=gg$group=="Group4"
  idx = idx1|idx2|idx3|idx4
  genes=as.character(gg$ENTREZID[idx])
  #====GO_KEGG_enrichment=====
  kegg1234=enrichment_analysis(geneList = genes, universe=universe, method = enrich_method,type = "KEGG",organism=organism,
                               pvalueCutoff = pvalue, plotTitle="KEGG: Group1&2&3&4", pAdjustMethod = adjust)
  bp1234=enrichment_analysis(geneList = genes, universe=universe, method = "ORT", type = "BP",organism=organism,
                             pvalueCutoff = pvalue, plotTitle="BP: Group1&2&3&4", pAdjustMethod = adjust)

  if(!is.null(filename)){
    if(!is.null(kegg1$enrichRes)){
      write.table(kegg1$enrichRes@result, file.path(out.dir,paste0("Group1_kegg_",filename,".txt")), sep="\t", row.names = F,col.names = T,quote=F)
      ggsave(kegg1$gridPlot,filename=file.path(out.dir,paste0("Group1_kegg_",filename,".png")),units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(bp1$enrichRes)){
      write.table(bp1$enrichRes@result, file.path(out.dir,paste0("Group1_bp_",filename,".txt")), sep="\t", row.names = F,col.names = T,quote=F)
      ggsave(bp1$gridPlot,filename=file.path(out.dir,paste0("Group1_bp_",filename,".png")),units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(kegg2$enrichRes)){
      write.table(kegg2$enrichRes@result, file.path(out.dir,paste0("Group2_kegg_",filename,".txt")), sep="\t", row.names = F,col.names = T,quote=F)
      ggsave(kegg2$gridPlot,filename=file.path(out.dir,paste0("Group2_kegg_",filename,".png")),units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(bp2$enrichRes)){
        write.table(bp2$enrichRes@result, file.path(out.dir,paste0("Group2_bp_",filename,".txt")), sep="\t", row.names = F,col.names = T,quote=F)
        ggsave(bp2$gridPlot,filename=file.path(out.dir,paste0("Group2_bp_",filename,".png")),units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(kegg3$enrichRes)){
      write.table(kegg3$enrichRes@result, file.path(out.dir,paste0("Group3_kegg_",filename,".txt")), sep="\t", row.names = F,col.names = T,quote=F)
      ggsave(kegg3$gridPlot,filename=file.path(out.dir,paste0("Group3_kegg_",filename,".png")),units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(bp3$enrichRes)){
      write.table(bp3$enrichRes@result, file.path(out.dir,paste0("Group3_bp_",filename,".txt")), sep="\t", row.names = F,col.names = T,quote=F)
      ggsave(bp3$gridPlot,filename=file.path(out.dir,paste0("Group3_bp_",filename,".png")),units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(kegg4$enrichRes)){
      write.table(kegg4$enrichRes@result, file.path(out.dir,paste0("Group4_kegg_",filename,".txt")), sep="\t", row.names = F,col.names = T,quote=F)
      ggsave(kegg4$gridPlot,filename=file.path(out.dir,paste0("Group4_kegg_",filename,".png")),units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(bp4$enrichRes)){
      write.table(bp4$enrichRes@result, file.path(out.dir,paste0("Group4_bp_",filename,".txt")), sep="\t", row.names = F,col.names = T,quote=F)
      ggsave(bp4$gridPlot,filename=file.path(out.dir,paste0("Group4_bp_",filename,".png")),units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(kegg13$enrichRes)){
      write.table(kegg13$enrichRes@result, file.path(out.dir,paste0("Group1&3_kegg_",filename,".txt")), sep="\t", row.names = F,col.names = T,quote=F)
      ggsave(kegg13$gridPlot,filename=file.path(out.dir,paste0("Group1&3_kegg_",filename,".png")),units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(bp13$enrichRes)){
      write.table(bp13$enrichRes@result, file.path(out.dir,paste0("Group1&3_bp_",filename,".txt")), sep="\t", row.names = F,col.names = T,quote=F)
      ggsave(bp13$gridPlot,filename=file.path(out.dir,paste0("Group1&3_bp_",filename,".png")),units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(kegg14$enrichRes)){
      write.table(kegg14$enrichRes@result, file.path(out.dir,paste0("Group1&4_kegg_",filename,".txt")), sep="\t", row.names = F,col.names = T,quote=F)
      ggsave(kegg14$gridPlot,filename=file.path(out.dir,paste0("Group1&4_kegg_",filename,".png")),units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(bp14$enrichRes)){
      write.table(bp14$enrichRes@result, file.path(out.dir,paste0("Group1&4_bp_",filename,".txt")), sep="\t", row.names = F,col.names = T,quote=F)
      ggsave(bp14$gridPlot,filename=file.path(out.dir,paste0("Group1&4_bp_",filename,".png")),units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(kegg23$enrichRes)){
      write.table(kegg23$enrichRes@result, file.path(out.dir,paste0("Group2&3_kegg_",filename,".txt")), sep="\t", row.names = F,col.names = T,quote=F)
      ggsave(kegg23$gridPlot,filename=file.path(out.dir,paste0("Group2&3_kegg_",filename,".png")),units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(bp23$enrichRes)){
      write.table(bp23$enrichRes@result, file.path(out.dir,paste0("Group2&3_bp_",filename,".txt")), sep="\t", row.names = F,col.names = T,quote=F)
      ggsave(bp23$gridPlot,filename=file.path(out.dir,paste0("Group2&3_bp_",filename,".png")),units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(kegg24$enrichRes)){
      write.table(kegg24$enrichRes@result, file.path(out.dir,paste0("Group2&4_kegg_",filename,".txt")), sep="\t", row.names = F,col.names = T,quote=F)
      ggsave(kegg24$gridPlot,filename=file.path(out.dir,paste0("Group2&4_kegg_",filename,".png")),units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(bp24$enrichRes)){
      write.table(bp24$enrichRes@result, file.path(out.dir,paste0("Group2&4_bp_",filename,".txt")), sep="\t", row.names = F,col.names = T,quote=F)
      ggsave(bp24$gridPlot,filename=file.path(out.dir,paste0("Group2&4_bp_",filename,".png")),units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(kegg1234$enrichRes)){
      write.table(kegg1234$enrichRes@result, file.path(out.dir,paste0("Group1&2&3&4_kegg_",filename,".txt")), sep="\t", row.names = F,col.names = T,quote=F)
      ggsave(kegg1234$gridPlot,filename=file.path(out.dir,paste0("Group1&2&3&4_kegg_",filename,".png")),units = "in",width=400/100,height =270/100 )
    }
    if(!is.null(bp1234$enrichRes)){
      write.table(bp1234$enrichRes@result, file.path(out.dir,paste0("Group1&2&3&4_bp_",filename,".txt")), sep="\t", row.names = F,col.names = T,quote=F)
      ggsave(bp1234$gridPlot,filename=file.path(out.dir,paste0("Group1&2&3&4_bp_",filename,".png")),units = "in",width=400/100,height =270/100 )
    }
  }
  return(list(kegg1=kegg1, kegg2=kegg2, kegg3=kegg3, kegg4=kegg4, kegg13=kegg13, kegg14=kegg14, kegg23=kegg23, kegg24=kegg24, kegg1234=kegg1234,
              bp1=bp1, bp2=bp2, bp3=bp3, bp4=bp4, bp13=bp13, bp14=bp14, bp23=bp23, bp24=bp24, bp1234=bp1234))

}

