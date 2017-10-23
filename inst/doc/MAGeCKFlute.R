## ----setup, echo=FALSE, results="hide"-----------------------------------
knitr::opts_chunk$set(tidy=FALSE, cache=TRUE,
                      dev="png",
                      message=FALSE, error=FALSE, warning=TRUE)

## ----quickStart1, eval=FALSE---------------------------------------------
#  ##Load MLE gene summary data
#  data("BRAF_mle.gene_summary")
#  ##Run the MAGeCKFlute MLE pipeline
#  FluteMLE(gene_summary=BRAF_mle.gene_summary,
#           ctrlname=c("D7_R1", "D7_R2"),
#           treatname=c("PLX7_R1","PLX7_R2"),
#           prefix="BRAF", organism="hsa")
#  
#  #All pipeline results are written into local directory
#  # "./BRAF_Flute_Results/", and all figures are integrated into file
#  # "BRAF_Flute.mle_summary.pdf".

## ----quickStart2, eval=FALSE---------------------------------------------
#  ##Load RRA gene summary
#  data("BRAF_rra.gene_summary")
#  ##Run the MAGeCKFlute RRA pipeline
#  FluteRRA(BRAF_rra.gene_summary, prefix="BRAF", organism="hsa")
#  
#  #All pipeline results are written into local directory
#  # "./BRAF_Flute_Results/" too, and all figures are integrated into file
#  # "BRAF_Flute.rra_summary.pdf".

## ----CheckMLERes---------------------------------------------------------
library(MAGeCKFlute)
data("BRAF_mle.gene_summary")
head(BRAF_mle.gene_summary)

## ----ReadBeta------------------------------------------------------------
gene_summary=BRAF_mle.gene_summary
ctrlName=c("D7_R1", "D7_R2")
treatName=c("PLX7_R1", "PLX7_R2")
#Read beta scores from gene summary matrix
dd=ReadBeta(gene_summary, ctrlName=ctrlName, treatName=treatName, organism="hsa")
head(dd)

## ----NormalizeBeta-------------------------------------------------------
dd_essential = NormalizeBeta(dd, method="cell_cycle")
head(dd_essential)

#OR
dd_loess = NormalizeBeta(dd, method="loess")
head(dd_loess)

## ----DistributeBeta------------------------------------------------------
Violin.plot(dd_essential, main="Cell cycle normalized")
Density.plot(dd_essential, main="Cell cycle normalized")
Density.diff(dd_essential, main="Cell cycle normalized")

#MAplot was used to evaluate the data quality of gene expression profile 
#in the past, here we can use it to evaluate the data quality of normalized
#beta score profile.
MA.plot(dd_essential, cex=1, main="Cell cycle normalized")

## ----EstimateCellCycle---------------------------------------------------
##Fitting beta score of all genes
CellCycleFit(dd_essential, ylab="Beta Score", main="Cell cycle normalized")

## ----selection-----------------------------------------------------------
ddAB_essential=GroupAB(dd_essential)
ScatterAB(ddAB_essential, main="Cell cycle normalized")
RankAB(ddAB_essential, main="Cell cycle normalized")

## ----EnrichAB, eval=FALSE------------------------------------------------
#  enrich_result = EnrichAB(ddAB_essential,pvalue=0.05,
#                           enrich_method = "Hypergeometric", organism="hsa")
#  print(enrich_result$keggA$gridPlot)
#  print(enrich_result$keggB$gridPlot)

## ----pathview, eval=FALSE,eval=FALSE-------------------------------------
#  genedata = dd_essential[,c("Control","Treatment")]
#  rownames(genedata)=dd_essential$ENTREZID
#  keggID = enrich_result$keggA$enrichRes@result$ID[1]
#  KeggPathwayView(gene.data = genedata, pathway.id = keggID, species="hsa")
#  #The pathway map will be located on current workspace

## ----9Square-------------------------------------------------------------
dd2=dd_essential[,c("Gene","Treatment","Control","ENTREZID")]
P2=Square.plot(dd2,main="Cell cycle normalized")
print(P2$p)

## ----EnrichSquare,eval=FALSE---------------------------------------------
#  enrich_result2 = EnrichSquare(P2$dd1,pvalue=0.05)
#  print(enrich_result2$kegg3$gridPlot)
#  print(enrich_result2$kegg4$gridPlot)

## ----pathview2, eval=FALSE,eval=FALSE------------------------------------
#  genedata = dd_essential[,c("Control","Treatment")]
#  rownames(genedata)=dd_essential$ENTREZID
#  keggID = enrich_result2$kegg1$enrichRes@result$ID[1]
#  KeggPathwayView(gene.data = genedata, pathway.id = keggID, species="hsa")

## ----CheckRRARes---------------------------------------------------------
data("BRAF_rra.gene_summary")
head(BRAF_rra.gene_summary)

## ----ReadRRA-------------------------------------------------------------
dd.rra = ReadRRA(BRAF_rra.gene_summary, organism="hsa")
head(dd.rra)

## ----selection2, eval=FALSE----------------------------------------------
#  ##Negative selection
#  universe=dd.rra$ENTREZID
#  genes = dd.rra[dd.rra$neg.fdr<0.05, "ENTREZID"]
#  kegg.neg=enrichment_analysis(geneList = genes, universe=universe,
#                               type = "KEGG", method="HyperGeometric",
#                               pvalueCutoff = 0.05, plotTitle="KEGG: neg")
#  bp.neg=enrichment_analysis(geneList = genes, universe=universe,
#                             method = "ORT", type = "BP",
#                             pvalueCutoff = 0.05, plotTitle="BP: neg")
#  print(kegg.neg$gridPlot)
#  print(bp.neg$gridPlot)
#  
#  ##Positive selection
#  universe=dd.rra$ENTREZID
#  genes = dd.rra[dd.rra$pos.fdr<0.05, "ENTREZID"]
#  kegg.pos=enrichment_analysis(geneList = genes, universe=universe,
#                               type = "KEGG", method="HyperGeometric",
#                               pvalueCutoff = 0.05, plotTitle="KEGG: pos")
#  bp.pos=enrichment_analysis(geneList = genes, universe=universe,
#                             method = "ORT", type = "BP",
#                             pvalueCutoff = 0.05, plotTitle="BP: pos")
#  print(kegg.pos$gridPlot)
#  print(bp.pos$gridPlot)
#  

## ----sessionInfo---------------------------------------------------------
sessionInfo()

