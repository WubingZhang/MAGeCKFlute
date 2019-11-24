## ----setup, echo=FALSE, fig.height=6, fig.width=9, dpi=300---------------
knitr::opts_chunk$set(tidy=FALSE, cache=TRUE,
                      dev="png", message=FALSE, error=FALSE, warning=TRUE)

## ----library, eval=TRUE--------------------------------------------------
library(MAGeCKFlute)
library(ggplot2)

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
BarView(countsummary, x = "Label", y = "GiniIndex",
        ylab = "Gini index", main = "Evenness of sgRNA reads")
countsummary$Missed = log10(countsummary$Zerocounts)
BarView(countsummary, x = "Label", y = "Missed", fill = "#394E80",
        ylab = "Log10 missed gRNAs", main = "Missed sgRNAs")
countsummary$Unmapped = countsummary$Reads - countsummary$Mapped
gg = data.table::melt(countsummary[, c("Label", "Mapped", "Unmapped")], id.vars = "Label")
gg$variable = factor(gg$variable, levels = c("Unmapped", "Mapped"))
p = BarView(gg, x = "Label", y = "value", fill = "variable", 
            position = "stack", xlab = NULL, ylab = "Reads", main = "Map ratio")
p + scale_fill_manual(values = c("#9BC7E9", "#1C6DAB"))


## ----CheckRRARes---------------------------------------------------------
library(MAGeCKFlute)

## ----ReadRRA-------------------------------------------------------------
dd.rra = ReadRRA(rra.gene_summary, score = "rra")
head(dd.rra)
dd.sgrna = ReadsgRRA(rra.sgrna_summary)

## ----selection1, fig.height=4, fig.width=7-------------------------------
dd.rra$LogFDR = -log10(dd.rra$FDR)
p1 = ScatterView(dd.rra, x = "LFC", y = "LogFDR", label = "id", model = "volcano", top = 5)
print(p1)

## ----rankrra, fig.height=4, fig.width=6----------------------------------
dd.rra$Rank = rank(dd.rra$LFC)
p2 = ScatterView(dd.rra, x = "LFC", y = "Rank", label = "id", 
                 top = 10, model = "rank", display_cut = T)
print(p2)

## ----scatter, fig.height=4, fig.width=6----------------------------------
dd.rra$RandomIndex = sample(1:nrow(dd.rra), nrow(dd.rra))
dd.rra = dd.rra[order(-dd.rra$LFC), ]
p2 = ScatterView(dd.rra, x = "RandomIndex", y = "LFC", label = "id",
                 y_cut = CutoffCalling(dd.rra$LFC,2), 
                 groups = "topright", toplabels = dd.rra$id[1:5])
p2 + ylim(0, NA)

## ----sgRNARank, fig.height=4, fig.width=7--------------------------------
p2 = sgRankView(dd.sgrna, top = 4, bottom = 4)
print(p2)

## ----enrich_rra----------------------------------------------------------
geneList= dd.rra$LFC
names(geneList) = dd.rra$id
enrich = EnrichAnalyzer(geneList = geneList, method = "GSEA", type = "GOBP")
EnrichedView(enrich)

## ----ReadBeta------------------------------------------------------------
data("mle.gene_summary")
dd=ReadBeta(mle.gene_summary)
head(dd)

## ----BatchRemove, fig.height=5, fig.width=10-----------------------------
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
ctrlname = "dmso"
treatname = "plx"
dd_essential = NormalizeBeta(dd, samples=c(ctrlname, treatname), method="cell_cycle")
head(dd_essential)

#OR
dd_loess = NormalizeBeta(dd, samples=c(ctrlname, treatname), method="loess")
head(dd_loess)

## ----DistributeBeta, fig.height=5, fig.width=8---------------------------
DensityView(dd_essential, samples=c(ctrlname, treatname))
ConsistencyView(dd_essential, ctrlname, treatname)

# Another option MAView
MAView(dd_essential, ctrlname, treatname)

## ----selection2, fig.height=5, fig.width=7-------------------------------
p1 = ScatterView(dd_essential, ctrlname, treatname, auto_cut_diag = TRUE)
print(p1)

## ----rank, fig.height=5, fig.width=7-------------------------------------
## Add column of 'diff'
dd_essential$Control = rowMeans(dd_essential[,ctrlname, drop = FALSE])
dd_essential$Treatment = rowMeans(dd_essential[,treatname, drop = FALSE])

dd_essential$Diff = dd_essential$Treatment - dd_essential$Control
dd_essential$Rank = rank(dd_essential$Diff)
p2 = ScatterView(dd_essential, x = "Diff", y = "Rank", label = "Gene", 
                 top = 10, model = "rank")
print(p2)

## ----Square, fig.height=7, fig.width=8-----------------------------------
p3 = ScatterView(dd_essential, x = "dmso", y = "plx", label = "Gene", 
                 model = "ninesquare", top = 5, display_cut = TRUE)
print(p3)

## ----EnrichSquare, fig.height=5, fig.width=9-----------------------------
#Get 9-square groups
Square9 = p3$data
idx=Square9$group=="topcenter"
geneList = Square9$Diff
names(geneList) = Square9$Gene[idx]
universe = Square9$Gene

# Enrichment analysis
kegg1 = EnrichAnalyzer(geneList = geneList, universe = universe)
EnrichedView(kegg1, top = 6, bottom = 0)

## ----pathview2, eval=FALSE-----------------------------------------------
#  genedata = dd_essential[, c("Control","Treatment")]
#  arrangePathview(genedata, pathways = "hsa01521", organism = "hsa", sub = NULL)

## ----sessionInfo---------------------------------------------------------
sessionInfo()

