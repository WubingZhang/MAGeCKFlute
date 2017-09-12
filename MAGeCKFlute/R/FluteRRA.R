#===read RRA results=====================================
FluteRRA <- function(gene_summary, prefix="Test", enrich_kegg="HyperGeometric", pvalueCutoff=0.05,
                     adjust="none", out.dir=paste0(prefix, "_Flute_Results"), workspace="."){
  #=========Prepare the running environment=========
  {
    loginfo("Create output dir and pdf file...")
    out.dir_sub=file.path(workspace, out.dir)
    dir.create(file.path(out.dir_sub), showWarnings=FALSE)
    dir.create(file.path(out.dir_sub,"RRA"), showWarnings=FALSE)

    output_pdf = paste0(prefix,"_Flute.rra_summary.pdf")
    pdf(file.path(out.dir_sub, output_pdf),width=9,height = 4)
  }

  #=========Input data=========
  loginfo("Read RRA result ...")
  dd = ReadRRA(gene_summary)

  #enrichment analysis
  {
    geneList = dd$neg.fdr
    names(geneList)=dd$ENTREZID
    genes = dd[dd$neg.fdr<pvalueCutoff, "ENTREZID"]

    kegg.neg=enrichment_analysis(geneList = geneList, genes=genes, method = enrich_kegg,type = "KEGG", pvalueCutoff = pvalueCutoff,
                               plotTitle="KEGG: neg",gridColour="#377eb8", pAdjustMethod = adjust)
    bp.neg=enrichment_analysis(geneList = geneList, genes=genes, method = "ORT", type = "BP", pvalueCutoff = pvalueCutoff,
                             plotTitle="BP: neg",gridColour="#377eb8", pAdjustMethod = adjust)

    grid.arrange(kegg.neg$gridPlot, bp.neg$gridPlot, ncol = 2)

    ggsave(kegg.neg$gridPlot,filename=file.path(out.dir_sub,"RRA/kegg.neg.png"),units = "in",width=400/100,height =270/100 )
    ggsave(bp.neg$gridPlot,filename=file.path(out.dir_sub,"RRA/bp.neg.png"),units = "in",width=400/100,height =270/100 )


    geneList = dd$pos.fdr
    names(geneList)=dd$ENTREZID
    genes = dd[dd$pos.fdr<pvalueCutoff, "ENTREZID"]

    kegg.pos=enrichment_analysis(geneList = geneList, genes=genes, method = enrich_kegg, type = "KEGG", pvalueCutoff = pvalueCutoff,
                             plotTitle="KEGG: pos",gridColour="#e41a1c", pAdjustMethod = adjust)
    bp.pos=enrichment_analysis(geneList = geneList, genes=genes, method = "ORT", type = "BP", pvalueCutoff = pvalueCutoff,
                           plotTitle="BP: pos",gridColour="#e41a1c", pAdjustMethod = adjust)
    # gse=enrichment_analysis(geneList = geneList, genes=genes, method = "GSEA", type = "KEGG", pvalueCutoff = pvalueCutoff,
    #                         plotTitle="GSEA: RRA",gridColour="#e41a1c", pAdjustMethod = adjust)
    grid.arrange(kegg.pos$gridPlot, bp.pos$gridPlot, ncol = 2)

    ggsave(kegg.pos$gridPlot,filename=file.path(out.dir_sub,"RRA/kegg.pos.png"),units = "in",width=400/100,height =270/100 )
    ggsave(bp.pos$gridPlot,filename=file.path(out.dir_sub,"RRA/bp.pos.png"),units = "in",width=400/100,height =270/100 )

  }
  dev.off()
}
