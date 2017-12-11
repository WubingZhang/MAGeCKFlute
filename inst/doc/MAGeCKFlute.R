## ----setup, echo=FALSE, results="hide"-----------------------------------
knitr::opts_chunk$set(tidy=FALSE, cache=TRUE,
                      dev="png",
                      message=FALSE, error=FALSE, warning=TRUE)

## ----quickStart1, eval=FALSE---------------------------------------------
#  ##Load gene summary data in MAGeCK MLE results
#  data("MLE_Data")
#  ##Run the downstream analysis pipeline for MAGeCK MLE
#  FluteMLE(MLE_Data, ctrlname=c("D7_R1", "D7_R2"), treatname=c("PLX7_R1","PLX7_R2"), prefix="BRAF", organism="hsa")
#  
#  #All pipeline results are written into local directory
#  # "./BRAF_Flute_Results/", and all figures are integrated into file "BRAF_Flute.mle_summary.pdf".

## ----quickStart2, eval=FALSE---------------------------------------------
#  ##Load gene summary data in MAGeCK RRA results
#  data("RRA_Data")
#  ##Run the downstream analysis pipeline for MAGeCK RRA
#  FluteRRA(RRA_Data, prefix="BRAF", organism="hsa")
#  
#  #All pipeline results are written into local directory
#  # "./BRAF_Flute_Results/" too, and all figures are integrated into file
#  # "BRAF_Flute.rra_summary.pdf".

## ----CheckMLERes---------------------------------------------------------
library(MAGeCKFlute)
data("MLE_Data")
head(MLE_Data)

## ----ReadBeta------------------------------------------------------------
gene_summary=MLE_Data
ctrlname = c("D7_R1", "D7_R2")
treatname = c("PLX7_R1", "PLX7_R2")
#Read beta scores from gene summary table in MAGeCK MLE results
dd=ReadBeta(gene_summary, organism="hsa")
head(dd)

## ----BatchRemove---------------------------------------------------------
##Before batch removal
HeatmapView(dd[,c(ctrlname, treatname)])
batchMat = data.frame(samples = c(ctrlname, treatname),
                      batch = c(1,2,1,2), cov = c(1,1,2,2))
dd1 = BatchRemove(dd, batchMat, cluster = FALSE)

##After batch removal
HeatmapView(dd[,c(ctrlname, treatname)])

## ----NormalizeBeta-------------------------------------------------------
dd_essential = NormalizeBeta(dd, method="cell_cycle")
head(dd_essential)

#OR
dd_loess = NormalizeBeta(dd, method="loess")
head(dd_loess)

## ----DistributeBeta------------------------------------------------------
ViolinView(dd_essential, samples=c(ctrlname, treatname), main="Cell cycle normalized")
DensityView(dd_essential, samples=c(ctrlname, treatname), main="Cell cycle normalized")
DensityDiffView(dd_essential, ctrlname, treatname, main="Cell cycle normalized")

#we can also use function 'MAView' to evaluate the data quality of normalized
#beta score profile.
MAView(dd_essential, ctrlname, treatname, cex=1, main="Cell cycle normalized")

## ----EstimateCellCycle---------------------------------------------------
##Fitting beta score of all genes
CellCycleView(dd_essential[,-2], ctrlname, ylab="Beta Score", main="Cell cycle normalized")

## ----selection-----------------------------------------------------------
p1 = ScatterView(dd_essential, ctrlname, treatname, main="Cell cycle normalized")
print(p1)

## ----rank----------------------------------------------------------------
## Add column of 'diff'
dd_essential$Control = rowMeans(dd_essential[,ctrlname])
dd_essential$Treatment = rowMeans(dd_essential[,treatname])
dd_essential$diff = dd_essential$Treatment - dd_essential$Control

p2 = RankView(dd_essential, main="Cell cycle normalized")
print(p2)

## ----EnrichAB, eval=TRUE-------------------------------------------------
## Get information of positive and negative selection genes
groupAB = p1$data
## select positive selection genes
idx1=groupAB$group=="up"
genes=as.character(groupAB$ENTREZID[idx1])
geneList=groupAB$diff[idx1]
names(geneList)=genes
universe=as.character(groupAB$ENTREZID)
## Do enrichment analysis using HGT method
keggA = enrichment_analysis(geneList = genes, universe = universe, method = "HGT",
                          type = "KEGG", organism = "human", plotTitle = "Positive selection")
#same as
kegg_A = enrich.HGT(genes, universe, type = "KEGG", organism = "human")
keggA_grid = EnrichedView(kegg_A@result, plotTitle = "Positive selection")

## look at the results
head(keggA$enrichRes@result)
print(keggA$gridPlot)
#should same as
head(kegg_A@result)
print(keggA_grid)


## ----GSEA----------------------------------------------------------------
## Do enrichment analysis using GSEA method
gseA = enrichment_analysis(geneList = geneList, method = "GSEA", type = "KEGG", 
                           organism="human", plotTitle="Positive selection")
#same as
gse_A = enrich.GSE(geneList, type = "KEGG", organism = "human")
gseA_grid = EnrichedGSEView(gse_A@result, plotTitle = "Positive selection")


head(gseA$enrichRes@result)
print(gseA$gridPlot)
#should same as
head(gse_A@result)
print(gseA_grid)

## ----pathview, eval=TRUE-------------------------------------------------
genedata = dd_essential[,c("Control","Treatment")]
rownames(genedata)=dd_essential$ENTREZID
keggID = keggA$enrichRes@result$ID[1]
#The pathway map will be located on current workspace
KeggPathwayView(gene.data = genedata, pathway.id = keggID, species="hsa")
##Read the figure into R
pngname=paste0(keggID, ".pathview.multi.png")
grid.arrange(grid::rasterGrob(png::readPNG(pngname)))
file.remove(pngname)

## ----9Square-------------------------------------------------------------
p3 = SquareView(dd_essential, main="Cell cycle normalized")
print(p3)

## ----EnrichSquare,eval=TRUE----------------------------------------------
##Get information of treatment-associated genes
Square9 = p3$data
##==select group1 genes in 9-Square
idx=Square9$group=="Group1"
genes=as.character(Square9$ENTREZID[idx])
universe=as.character(Square9$ENTREZID)
#====KEGG_enrichment=====
kegg1=enrichment_analysis(geneList = genes, universe = universe, 
                          type = "KEGG",plotTitle = "KEGG: Group1")
## look at the results
head(kegg1$enrichRes@result)
print(kegg1$gridPlot)

## ----pathview2, eval=FALSE,eval=TRUE-------------------------------------
genedata = dd_essential[, c("Control","Treatment")]
rownames(genedata) = dd_essential$ENTREZID
keggID = kegg1$enrichRes@result$ID[1]
KeggPathwayView(gene.data = genedata, pathway.id = keggID, species="hsa")
##Read the figure into R
pngname=paste0(keggID, ".pathview.multi.png")
grid.arrange(grid::rasterGrob(png::readPNG(pngname)))
file.remove(pngname)

## ----CheckRRARes---------------------------------------------------------
data("RRA_Data")
head(RRA_Data)

## ----ReadRRA-------------------------------------------------------------
dd.rra = ReadRRA(RRA_Data, organism="hsa")
head(dd.rra)

## ----selection2, eval=FALSE----------------------------------------------
#  ##Negative selection
#  universe=dd.rra$ENTREZID
#  genes = dd.rra[dd.rra$neg.fdr<0.05, "ENTREZID"]
#  kegg.neg=enrichment_analysis(geneList = genes, universe=universe,
#                               type = "KEGG", plotTitle="KEGG: neg")
#  bp.neg=enrichment_analysis(geneList = genes, universe=universe,
#                             type = "BP", plotTitle="BP: neg")
#  print(kegg.neg$gridPlot)
#  print(bp.neg$gridPlot)
#  
#  ##Positive selection
#  universe=dd.rra$ENTREZID
#  genes = dd.rra[dd.rra$pos.fdr<0.05, "ENTREZID"]
#  kegg.pos=enrichment_analysis(geneList = genes, universe=universe,
#                               type = "KEGG", plotTitle="KEGG: pos")
#  bp.pos=enrichment_analysis(geneList = genes, universe=universe,
#                             type = "BP", plotTitle="BP: pos")
#  print(kegg.pos$gridPlot)
#  print(bp.pos$gridPlot)
#  

## ----sessionInfo---------------------------------------------------------
sessionInfo()

