## ----setup, echo=FALSE, fig.height=4, fig.width=20, dpi=150--------------
knitr::opts_chunk$set(tidy=FALSE, cache=TRUE,
                      dev="png", message=FALSE, error=FALSE, warning=TRUE)

## ----load----------------------------------------------------------------
library(MAGeCKFlute)
data(rra.gene_summary)
df = ReadRRA(rra.gene_summary)
genelist= df$Score
names(genelist) = df$id

## ----HGT-----------------------------------------------------------------
# Alternative functions EnrichAnalyzer and enrich.HGT.
hgtRes1 = EnrichAnalyzer(genelist[genelist< -1], method = "HGT")
head(hgtRes1@result)
hgtRes2 = enrich.HGT(genelist[genelist< -1])
head(hgtRes2@result)

## ----ORT-----------------------------------------------------------------
# Alternative functions EnrichAnalyzer and enrich.ORT.
ortRes1 = EnrichAnalyzer(genelist[genelist< -1], method = "ORT")
head(ortRes1@result)
ortRes2 = enrich.ORT(genelist[genelist< -1])
head(ortRes2@result)

## ----GSE-----------------------------------------------------------------
# Alternative functions EnrichAnalyzer and enrich.GSE.
gseRes1 = EnrichAnalyzer(genelist, method = "GSEA")
head(gseRes1@result)
gseRes2 = enrich.GSE(genelist)
head(gseRes2@result)

## ----barview, fig.height=4, fig.width=9, dpi=150-------------------------
require(ggplot2)
df = hgtRes1@result
df$logFDR = -log10(df$p.adjust)
p = BarView(df[1:5,], "Description", 'logFDR')
p = p + labs(x = NULL) + coord_flip()
p

# Or use function barplot from enrichplot package
barplot(hgtRes1, showCategory = 5)

## ---- fig.height=4, fig.width=9, dpi=150---------------------------------
EnrichedView(hgtRes1, bottom = 5, mode = 1)
EnrichedView(hgtRes1, bottom = 5, mode = 2)
dotplot(hgtRes1, showCategory = 5)

## ---- fig.height=7, fig.width=9, dpi=150---------------------------------
hgtRes1@result$geneID = hgtRes1@result$geneName
#cnetplot
cnetplot(hgtRes1, 2)

## ---- fig.height=5, fig.width=12, dpi=150--------------------------------
heatplot(hgtRes1, showCategory = 3, foldChange=genelist)
emapplot(hgtRes1, layout="kk")

## ---- fig.height=4, fig.width=9, dpi=150---------------------------------
#gseaplot
gseaplot(gseRes1, geneSetID = 1, title = gseRes1$Description[1])
gseaplot(gseRes1, geneSetID = 1, by = "runningScore", title = gseRes1$Description[1])
gseaplot(gseRes1, geneSetID = 1, by = "preranked", title = gseRes1$Description[1])
#or
gseaplot2(gseRes1, geneSetID = 1:3)

## ----pathway, fig.height=4, fig.width=7, dpi=150-------------------------
## KEGG and REACTOME pathways
enrich = EnrichAnalyzer(geneList = genelist[genelist< -1], type = "KEGG+REACTOME")
EnrichedView(enrich, bottom = 5)
## Only KEGG pathways
enrich = EnrichAnalyzer(geneList = genelist[genelist< -1], type = "KEGG")
EnrichedView(enrich, bottom = 5)
## Gene ontology
enrichGo = EnrichAnalyzer(genelist[genelist< -1], type = "GOBP+GOMF")
EnrichedView(enrichGo, bottom = 5)

## ---- fig.height=4, fig.width=7, dpi=150---------------------------------
enrichPro = EnrichAnalyzer(genelist[genelist< -1], type = "CORUM+CPX")
EnrichedView(enrichPro, bottom = 5)

## ---- fig.height=4, fig.width=7, dpi=150---------------------------------
enrichComb = EnrichAnalyzer(genelist[genelist< -1], type = "GOBP+KEGG")
EnrichedView(enrichComb, bottom = 5)

## ----limit, fig.height=4, fig.width=7, dpi=150---------------------------
enrich = EnrichAnalyzer(genelist[genelist< -1], type = "GOBP", limit = c(1, 80))
EnrichedView(enrich, bottom = 5)

## ----filter1, fig.height=4, fig.width=7, dpi=150-------------------------
enrich1 = EnrichAnalyzer(genelist[genelist< -1], type = "GOMF+GOBP")
enrich2 = EnrichAnalyzer(genelist[genelist< -1], type = "GOMF+GOBP", filter = TRUE)
enrich3 = EnrichedFilter(enrich1)
EnrichedView(enrich1, bottom = 15)
EnrichedView(enrich2, bottom = 15)
EnrichedView(enrich3, bottom = 15)

## ----sessionInfo---------------------------------------------------------
sessionInfo()

