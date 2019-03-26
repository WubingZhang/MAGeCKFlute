## ----setup, echo=FALSE, fig.height=6, fig.width=9, dpi=300---------------
knitr::opts_chunk$set(tidy=FALSE, cache=TRUE,
                      dev="png", message=FALSE, error=FALSE, warning=TRUE)

## ----library, eval=TRUE--------------------------------------------------
library(MAGeCKFlute)

## ----quickStart2, eval=FALSE---------------------------------------------
#  ##Load gene summary data in MAGeCK RRA results
#  data("rra.gene_summary")
#  data("rra.sgrna_summary")
#  ##Run the downstream analysis pipeline for MAGeCK RRA
#  FluteRRA(rra.gene_summary, rra.sgrna_summary, prefix="RRA", organism="hsa")

## ----quickStart1, eval=FALSE---------------------------------------------
#  ## Load gene summary data in MAGeCK MLE results
#  data("mle.gene_summary")
#  ## Run the downstream analysis pipeline for MAGeCK MLE
#  FluteMLE(mle.gene_summary, ctrlname="dmso", treatname="plx", prefix="MLE", organism="hsa")

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

## ----CheckRRARes---------------------------------------------------------
library(MAGeCKFlute)
data("rra.gene_summary")
head(rra.gene_summary)

## ----ReadRRA-------------------------------------------------------------
dd.rra = ReadRRA(rra.gene_summary, organism = "hsa")
head(dd.rra)
dd.sgrna = ReadsgRRA(rra.sgrna_summary)

## ----selection1, fig.height=4, fig.width=7-------------------------------
p1 = VolcanoView(dd.rra, x = "LFC", y = "FDR", Label = "Official")
print(p1)

## ----rankrra, fig.height=4, fig.width=6----------------------------------
geneList= dd.rra$LFC
names(geneList) = dd.rra$Official
p2 = RankView(geneList, top = 10, bottom = 10)
print(p2)

## ----sgRNARank, fig.height=4, fig.width=7--------------------------------
p2 = sgRankView(dd.sgrna, top = 0, bottom = 0, gene = levels(p1$data$Label))
print(p2)

## ----enrich_rra----------------------------------------------------------
universe = dd.rra$EntrezID
geneList= dd.rra$LFC
names(geneList) = universe

enrich = EnrichAnalyzer(geneList = geneList, method = "GSEA", type = "GOMF+GOCC+GOBP", limit = c(2, 100))

## ----enrichedGeneView, fig.height=5, fig.width=15------------------------
EnrichedGeneView(slot(enrich, "result"), geneList, keytype = "Entrez")
EnrichedView(slot(enrich, "result"))

## ------------------------------------------------------------------------
enrich = EnrichAnalyzer(geneList = geneList, method = "GSEA", type = "GOMF+GOCC+GOBP", limit = c(2, 100), filter = FALSE)
enrich2 = EnrichedFilter(enrich)
EnrichedView(enrich2)

## ----CheckMLERes---------------------------------------------------------
library(MAGeCKFlute)
data("mle.gene_summary")
head(mle.gene_summary)

## ----ReadBeta------------------------------------------------------------
data("mle.gene_summary")
ctrlname = c("dmso")
treatname = c("plx")
#Read beta scores from gene summary table in MAGeCK MLE results
dd=ReadBeta(mle.gene_summary, organism="hsa")
head(dd)

## ----BatchRemove, fig.height=5, fig.width=6------------------------------
##Before batch removal
edata = matrix(c(rnorm(2000, 5), rnorm(2000, 8)), 1000)
colnames(edata) = paste0("s", 1:4)
HeatmapView(cor(edata))

## After batch removal
batchMat = data.frame(sample = colnames(edata), batch = rep(1:2, each = 2))
edata1 = BatchRemove(edata, batchMat)
head(edata1$data)
print(edata1$p)

## ----NormalizeBeta-------------------------------------------------------
dd_essential = NormalizeBeta(dd, samples=c(ctrlname, treatname), method="cell_cycle")
head(dd_essential)

#OR
dd_loess = NormalizeBeta(dd, samples=c(ctrlname, treatname), method="loess")
head(dd_loess)

## ----DistributeBeta, fig.height=5, fig.width=8---------------------------
ViolinView(dd_essential, samples=c(ctrlname, treatname))
DensityView(dd_essential, samples=c(ctrlname, treatname))
DensityDiffView(dd_essential, ctrlname, treatname)

#we can also use the function 'MAView' to evaluate the data quality of normalized
#beta score profile.
MAView(dd_essential, ctrlname, treatname)

## ----EstimateCellCycle, fig.height=5, fig.width=8------------------------
##Fitting beta score of all genes
CellCycleView(dd_essential, ctrlname, treatname)

## ----selection2, fig.height=5, fig.width=7-------------------------------
p1 = ScatterView(dd_essential, ctrlname, treatname)
print(p1)

## ----rank, fig.height=5, fig.width=7-------------------------------------
## Add column of 'diff'
dd_essential$Control = rowMeans(dd_essential[,ctrlname, drop = FALSE])
dd_essential$Treatment = rowMeans(dd_essential[,treatname, drop = FALSE])

rankdata = dd_essential$Treatment - dd_essential$Control
names(rankdata) = dd_essential$Gene
p2 = RankView(rankdata)
print(p2)

## ----EnrichAB, fig.height=5, fig.width=10--------------------------------
## Get information of positive and negative selection genes
groupAB = p1$data
geneList = groupAB$diff; names(geneList) = rownames(groupAB)
## Do enrichment analysis for positive selection genes.
idx1 = groupAB$group=="up"
hgtA = EnrichAnalyzer(geneList[idx1], method = "HGT", universe = rownames(groupAB))
hgtA_grid = EnrichedView(slot(hgtA, "result"))

## look at the results
head(slot(hgtA, "result"))
print(hgtA_grid)

## ----GSEA, fig.height=5, fig.width=10------------------------------------
## Do enrichment analysis using GSEA method
gseA = EnrichAnalyzer(geneList, method = "GSEA", type = "KEGG", limit = c(2, 150))
gseA_grid = EnrichedView(gseA)
print(gseA_grid)

## ----pathview, fig.height=10, fig.width=20-------------------------------
genedata = dd_essential[,c("Control","Treatment")]
keggID = gsub("KEGG_", "", slot(gseA, "result")$ID[1])
arrangePathview(genedata, pathways = keggID, organism = "hsa", sub = NULL)

## ----Square, fig.height=7, fig.width=8-----------------------------------
p3 = SquareView(dd_essential, label = "Gene")
print(p3)

## ----EnrichSquare, fig.height=5, fig.width=9-----------------------------
#Get 9-square groups
Square9 = p3$data
idx=Square9$group=="topcenter"
geneList = (Square9$y - Square9$x)[idx]
names(geneList) = rownames(Square9)[idx]
universe=rownames(Square9)
# Enrichment analysis
kegg1 = EnrichAnalyzer(geneList = geneList, universe = universe)
EnrichedView(kegg1, top = 10, bottom = 0)

## ----pathview2, eval=FALSE-----------------------------------------------
#  genedata = dd_essential[, c("Control","Treatment")]
#  arrangePathview(genedata, pathways = "hsa01521", organism = "hsa", sub = NULL)

## ----sessionInfo---------------------------------------------------------
sessionInfo()

