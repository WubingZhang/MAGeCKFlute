## ----setup, echo=FALSE, fig.height=6, fig.width=9, dpi=300---------------
knitr::opts_chunk$set(tidy=FALSE, cache=TRUE,
                      dev="png", message=FALSE, error=FALSE, warning=TRUE)

## ----library, eval=TRUE--------------------------------------------------
library(MAGeCKFlute)
library(ggplot2)

## ----quickStart, eval=FALSE----------------------------------------------
#  ##Load gene summary data in MAGeCK RRA results
#  data("rra.gene_summary")
#  data("rra.sgrna_summary")
#  ##Run the downstream analysis pipeline for MAGeCK RRA
#  FluteRRA(rra.gene_summary, rra.sgrna_summary, proj="PLX", organism="hsa")

## ----quickStart1, eval=FALSE---------------------------------------------
#  data("mle.gene_summary")
#  ## Run the downstream analysis pipeline for MAGeCK MLE
#  FluteMLE(mle.gene_summary, treatname="plx", ctrlname="dmso", proj="PLX", organism="hsa")

## ----quickStart2, eval=FALSE---------------------------------------------
#  ## Take Depmap screen as control
#  FluteMLE(mle.gene_summary, treatname="plx", ctrlname="Depmap", proj="PLX", organism="hsa", incorporateDepmap = TRUE)
#  ## If a specific cell line is preferred, you can look at the similarity of your screen with Depmap screens by using function ResembleDepmap.
#  depmap_similarity = ResembleDepmap(mle.gene_summary, symbol = "Gene", score = "plx.beta")
#  FluteMLE(mle.gene_summary, treatname="plx", ctrlname="Depmap", proj="PLX", organism="hsa", incorporateDepmap = TRUE, cell_lines = rownames(depmap_similarity)[1], lineages = "All")

## ----CheckCountSummary---------------------------------------------------
data("countsummary")
head(countsummary)

## ----CountQC, fig.height=5, fig.width=7----------------------------------
# Gini index
BarView(countsummary, x = "Label", y = "GiniIndex",
        ylab = "Gini index", main = "Evenness of sgRNA reads")

# Missed sgRNAs
countsummary$Missed = log10(countsummary$Zerocounts)
BarView(countsummary, x = "Label", y = "Missed", fill = "#394E80",
        ylab = "Log10 missed gRNAs", main = "Missed sgRNAs")

# Read mapping
MapRatesView(countsummary)
# Or
countsummary$Unmapped = countsummary$Reads - countsummary$Mapped
gg = data.table::melt(countsummary[, c("Label", "Mapped", "Unmapped")], id.vars = "Label")
gg$variable = factor(gg$variable, levels = c("Unmapped", "Mapped"))
p = BarView(gg, x = "Label", y = "value", fill = "variable", 
            position = "stack", xlab = NULL, ylab = "Reads", main = "Map ratio")
p + scale_fill_manual(values = c("#9BC7E9", "#1C6DAB"))


## ----CheckRRARes---------------------------------------------------------
library(MAGeCKFlute)
data("rra.gene_summary")
head(rra.gene_summary)

## ----ReadRRA-------------------------------------------------------------
dd.rra = ReadRRA(rra.gene_summary)
head(dd.rra)
dd.sgrna = ReadsgRRA(rra.sgrna_summary)

## ------------------------------------------------------------------------
depmap_similarity = ResembleDepmap(dd.rra, symbol = "id", score = "Score", lineages = "All", method = "pearson")
head(depmap_similarity)

## ----omitessential-------------------------------------------------------
dd.rra = OmitCommonEssential(dd.rra)
dd.sgrna = OmitCommonEssential(dd.sgrna, symbol = "Gene")
# Compute the similarity with Depmap screens based on subset genes
depmap_similarity = ResembleDepmap(dd.rra, symbol = "id", score = "Score", lineages = "All", method = "pearson")
head(depmap_similarity)

## ----selection1, fig.height=4, fig.width=7-------------------------------
dd.rra$LogFDR = -log10(dd.rra$FDR)
p1 = ScatterView(dd.rra, x = "Score", y = "LogFDR", label = "id", model = "volcano", top = 5)
print(p1)

# Or
p2 = VolcanoView(dd.rra, x = "Score", y = "FDR", Label = "id")
print(p2)


## ----rankrra, fig.height=6, fig.width=4----------------------------------
dd.rra$Rank = rank(dd.rra$Score)
p1 = ScatterView(dd.rra, x = "Rank", y = "Score", label = "id", 
                 top = 5, auto_cut_y = TRUE, groups = c("top", "bottom"))
print(p1)
# Label interested hits using parameter `toplabels` (in ScatterView) and `genelist` (in RankView).
ScatterView(dd.rra, x = "Rank", y = "Score", label = "id", auto_cut_y = TRUE, 
            groups = c("top", "bottom"), toplabels = c("NF1", "EP300", "NF2"))

## ----rankrra2, fig.height=5, fig.width=6---------------------------------
# Or
geneList= dd.rra$Score
names(geneList) = dd.rra$id
p2 = RankView(geneList, top = 5, bottom = 10) + ylab("Score")
print(p2)
RankView(geneList, top = 0, bottom = 0, genelist = c("NF1", "EP300", "NF2")) + ylab("Score")

## ----scatter, fig.height=4, fig.width=6----------------------------------
dd.rra$RandomIndex = sample(1:nrow(dd.rra), nrow(dd.rra))
dd.rra = dd.rra[order(-dd.rra$Score), ]
p2 = ScatterView(dd.rra, x = "RandomIndex", y = "Score", label = "id",
                 y_cut = CutoffCalling(dd.rra$Score,2), 
                 groups = "top", toplabels = dd.rra$id[1:5])
p2 + ylim(0, NA)

## ----sgRNARank, fig.height=4, fig.width=7--------------------------------
p2 = sgRankView(dd.sgrna, top = 4, bottom = 4)
print(p2)

## ----enrich_rra, fig.height=4, fig.width=9-------------------------------
geneList= dd.rra$Score
names(geneList) = dd.rra$id
enrich = EnrichAnalyzer(geneList = geneList, method = "GSEA", type = "GOBP")
EnrichedView(enrich)

## ----ReadBeta------------------------------------------------------------
data("mle.gene_summary")
dd=ReadBeta(mle.gene_summary)
head(dd)

## ----BatchRemove, fig.height=6, fig.width=9------------------------------
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
dd_essential$Control = rowMeans(dd_essential[,ctrlname, drop = FALSE])
dd_essential$Treatment = rowMeans(dd_essential[,treatname, drop = FALSE])

p1 = ScatterView(dd_essential, "Control", "Treatment", groups = c("top", "bottom"), auto_cut_diag = TRUE, display_cut = TRUE, toplabels = c("NF1", "NF2", "EP300"))
print(p1)

## ----rank, fig.height=5, fig.width=7-------------------------------------
dd_essential$Diff = dd_essential$Treatment - dd_essential$Control
dd_essential$Rank = rank(dd_essential$Diff)
p1 = ScatterView(dd_essential, x = "Diff", y = "Rank", label = "Gene", 
                 top = 5, model = "rank")
print(p1)

# Or
rankdata = dd_essential$Treatment - dd_essential$Control
names(rankdata) = dd_essential$Gene
RankView(rankdata)

## ----Square, fig.height=6, fig.width=8-----------------------------------
p1 = ScatterView(dd_essential, x = "dmso", y = "plx", label = "Gene", 
                 model = "ninesquare", top = 5, display_cut = TRUE, force = 2)
print(p1)

# Or
p2 = SquareView(dd_essential, label = "Gene", 
                x_cutoff = CutoffCalling(dd_essential$Control, 2), 
                y_cutoff = CutoffCalling(dd_essential$Treatment, 2))
print(p2)

## ----EnrichSquare, fig.height=4, fig.width=9-----------------------------
# 9-square groups
Square9 = p1$data
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

