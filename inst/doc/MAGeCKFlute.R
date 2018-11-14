## ----setup, echo=FALSE, fig.height=6, fig.width=9, dpi=300---------------
knitr::opts_chunk$set(tidy=FALSE, cache=TRUE,
                      dev="png", message=FALSE, error=FALSE, warning=TRUE)

## ----quickStart1, eval=FALSE---------------------------------------------
#  library(MAGeCKFlute)
#  ##Load gene summary data in MAGeCK MLE results
#  data("MLE_Data")
#  ##Run the downstream analysis pipeline for MAGeCK MLE
#  FluteMLE(MLE_Data, ctrlname=c("D7_R1", "D7_R2"), treatname=c("PLX7_R1","PLX7_R2"), prefix="BRAF_", organism="hsa")

## ----quickStart2, eval=FALSE---------------------------------------------
#  ##Load gene summary data in MAGeCK RRA results
#  data("RRA_Data")
#  ##Run the downstream analysis pipeline for MAGeCK RRA
#  FluteRRA(RRA_Data, prefix="BRAF", organism="hsa")

## ----CheckRRARes---------------------------------------------------------
library(MAGeCKFlute)
data("RRA_Data")
head(RRA_Data)

## ----ReadRRA-------------------------------------------------------------
dd.rra = ReadRRA(RRA_Data, organism = "hsa")
head(dd.rra)

## ----selection1, fig.height=4, fig.width=7-------------------------------
RRA_Data$fdr = rowMin(as.matrix(RRA_Data[, c("neg.fdr", "pos.fdr")]))
p1 = VolcanoView(RRA_Data, x = "neg.lfc", y = "fdr", Label = "id")
print(p1)

## ----neg3, fig.height=5, fig.width=9-------------------------------------
universe=dd.rra$ENTREZID
idx = dd.rra$neg.fdr<0.05
geneList= -log10(dd.rra[idx, "neg.fdr"])
names(geneList) = dd.rra$ENTREZID[idx]

kegg.neg = enrichment_analysis(geneList = geneList, universe=universe, 
                             type = "KEGG", plotTitle="KEGG: neg")
bp.neg = enrichment_analysis(geneList = geneList, universe=universe, 
                           type = "GOBP", plotTitle="GOBP: neg")
print(kegg.neg$gridPlot)
print(bp.neg$gridPlot)

## ----pos3, fig.height=5, fig.width=10------------------------------------
universe = dd.rra$ENTREZID
idx = dd.rra$pos.fdr<0.05
geneList = -log10(dd.rra[idx, "pos.fdr"])
names(geneList) = dd.rra$ENTREZID[idx]
kegg.pos=enrichment_analysis(geneList = geneList, universe=universe, 
                             type = "KEGG", plotTitle="KEGG: pos")
bp.pos=enrichment_analysis(geneList = geneList, universe=universe, 
                           type = "GOBP", plotTitle="GOBP: pos")
print(kegg.pos$gridPlot)
print(bp.pos$gridPlot)

## ----CheckCountSummary---------------------------------------------------
data("countsummary")
head(countsummary)

## ----CountQC, fig.height=5, fig.width=7----------------------------------
MapRatesView(countsummary)
IdentBarView(countsummary, x = "Label", y = "GiniIndex", 
             ylab = "Gini index", main = "Evenness of sgRNA reads")
countsummary$Missed = log10(countsummary$Zerocounts)
IdentBarView(countsummary, x = "Label", y = "Missed", fill = "#394E80",
             ylab = "Log10 missed gRNAs", main = "Missed sgRNAs")

## ----CheckMLERes---------------------------------------------------------
library(MAGeCKFlute)
data("MLE_Data")
head(MLE_Data)

## ----ReadBeta------------------------------------------------------------
data("MLE_Data")
gene_summary = MLE_Data
ctrlname = c("D7_R1", "D7_R2")
treatname = c("PLX7_R1", "PLX7_R2")
#Read beta scores from gene summary table in MAGeCK MLE results
dd=ReadBeta(gene_summary, organism="hsa")
head(dd)

## ----BatchRemove, fig.height=5, fig.width=6------------------------------
##Before batch removal
HeatmapView(dd[,c(ctrlname, treatname)])
batchMat = data.frame(samples = c(ctrlname, treatname), batch = c(1,2,1,2), cov = c(1,1,2,2))
dd1 = BatchRemove(dd[,c(ctrlname, treatname)], batchMat)$data

## After batch removal
HeatmapView(dd1[,c(ctrlname, treatname)])

## ----NormalizeBeta-------------------------------------------------------
dd_essential = NormalizeBeta(dd, samples=c(ctrlname, treatname), method="cell_cycle")
head(dd_essential)

#OR
dd_loess = NormalizeBeta(dd, samples=c(ctrlname, treatname), method="loess")
head(dd_loess)

## ----DistributeBeta, fig.height=5, fig.width=8---------------------------
ViolinView(dd_essential, samples=c(ctrlname, treatname), main="Cell cycle normalized")
DensityView(dd_essential, samples=c(ctrlname, treatname), main="Cell cycle normalized")
DensityDiffView(dd_essential, ctrlname, treatname, main="Cell cycle normalized")

#we can also use the function 'MAView' to evaluate the data quality of normalized
#beta score profile.
MAView(dd_essential, ctrlname, treatname, cex=1, main="Cell cycle normalized")

## ----EstimateCellCycle, fig.height=5, fig.width=8------------------------
##Fitting beta score of all genes
CellCycleView(dd_essential, ctrlname, treatname, main="Cell cycle normalized")

## ----selection2, fig.height=5, fig.width=7-------------------------------
p1 = ScatterView(dd_essential, ctrlname, treatname, main="Cell cycle normalized")
print(p1)

## ----rank, fig.height=5, fig.width=7-------------------------------------
## Add column of 'diff'
dd_essential$Control = rowMeans(dd_essential[,ctrlname])
dd_essential$Treatment = rowMeans(dd_essential[,treatname])

rankdata = dd_essential$Treatment - dd_essential$Control
names(rankdata) = dd_essential$Gene
p2 = RankView(rankdata, main="Cell cycle normalized")
print(p2)

## ----EnrichAB, fig.height=5, fig.width=10--------------------------------
## Get information of positive and negative selection genes
groupAB = p1$data
## select positive selection genes
idx1=groupAB$group=="up"
genes=rownames(groupAB)[idx1]
geneList=groupAB$diff[idx1]
names(geneList)=genes
geneList = sort(geneList, decreasing = TRUE)
universe=rownames(groupAB)
## Do enrichment analysis using HGT method
keggA = enrich.HGT(geneList[1:100], universe, organism = "hsa", limit = c(3, 50))
keggA_grid = EnrichedGSEView(keggA@result, plotTitle = "Positive selection")

## look at the results
head(keggA@result)
print(keggA_grid)


## ----GSEA, fig.height=5, fig.width=10------------------------------------
## Do enrichment analysis using GSEA method
gseA = enrich.GSE(geneList, type = "KEGG", organism = "hsa", pvalueCutoff = 1)
gseA_grid = EnrichedGSEView(gseA@result, plotTitle = "Positive selection")

#should same as
head(gseA@result)
print(gseA_grid)

## ----pathview, fig.height=10, fig.width=20-------------------------------
genedata = dd_essential[,c("Control","Treatment")]
keggID = gsub("KEGG_", "", gseA@result$ID[1])
#The pathway map will be located on current workspace
KeggPathwayView(gene.data = genedata, pathway.id = keggID, species = "hsa")
##Read the figure into R
pngname=paste0(keggID, ".pathview.multi.png")
grid.arrange(grid::rasterGrob(png::readPNG(pngname)))
file.remove(paste0(keggID, c(".pathview.multi.png", ".png", ".xml")))

## ----Square, fig.height=7, fig.width=8-----------------------------------
p3 = SquareView(dd_essential, label = "Gene", main="Cell cycle normalized")
print(p3)

## ----EnrichSquare, fig.height=5, fig.width=9-----------------------------
##Get information of treatment-associated genes
Square9 = p3$data
##==select group1 genes in 9-Square
idx=Square9$group=="Group1"
geneList = (Square9$Treatment - Square9$Control)[idx]
names(geneList) = rownames(Square9)[idx]
universe=rownames(Square9)
#====KEGG_enrichment=====
kegg1=enrich.ORT(geneList = geneList, universe = universe, type = "KEGG", limit = c(3, 50))
## look at the results
head(kegg1@result)
EnrichedGSEView(kegg1@result)

## ----pathview2, eval=FALSE-----------------------------------------------
#  genedata = dd_essential[, c("Control","Treatment")]
#  keggID = gsub("KEGG_", "", kegg1@result$ID[1])
#  KeggPathwayView(gene.data = genedata, pathway.id = keggID, species="hsa")
#  ##Read the figure into R
#  pngname=paste0(keggID, ".pathview.multi.png")
#  grid.arrange(grid::rasterGrob(png::readPNG(pngname)))
#  file.remove(paste0(keggID, c(".pathview.multi.png", ".png", ".xml")))

## ----sessionInfo---------------------------------------------------------
sessionInfo()

