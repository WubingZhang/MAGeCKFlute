#' Visualize selected genes in enriched genesets
#'
#' @docType methods
#' @name EnrichedGeneView
#' @rdname EnrichedGeneView
#'
#' @param enrichment A data frame of enrichment result.
#' @param geneList The geneList used in enrichment analysis.
#' @param keytype "Entrez" or "Symbol".
#' @param gene_cutoff A tow-length numeric vector, specifying cutoff for negative and positive selections.
#' @param top An integer, specifying the number of top enriched terms to show.
#' @param bottom An integer, specifying the number of bottom enriched terms to show.
#' @param charLength Integer, specifying max length of enriched term name to show as coordinate lab.
#' @param filename Figure file name to create on disk. Default filename="NULL", which means no output.
#' @param width As in ggsave.
#' @param height As in ggsave.
#' @param ... Other available parameters in ggsave.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#' @examples
#' data(geneList, package = "DOSE")
#' enrichRes <- enrich.GSE(geneList)
#' EnrichedGeneView(enrichment=as.data.frame(enrichRes), geneList, keytype = "Entrez")
#' @export

EnrichedGeneView=function(enrichment, geneList, keytype = "Symbol", gene_cutoff = c(-log2(1.5), log2(1.5)),
                          top = 5, bottom = 5, charLength = 40, filename = NULL, width = 7, height = 5, ...){
  # No enriched pathways
  if(is.null(enrichment) || nrow(enrichment)==0){
    p1 = noEnrichPlot("No enriched terms")
    if(!is.null(filename)){
      ggsave(plot=p1,filename=filename, units = "in", width=width, height=height, ...)
    }
    return(p1)
  }
  # Reorder enrichment table
  enrichment$p.adjust = p.adjust(enrichment$pvalue, "BH")
  enrichment$logP = round(-log10(enrichment$p.adjust), 1)
  enrichment = enrichment[!is.na(enrichment$ID),]
  enrichment=enrichment[!duplicated(enrichment$Description),]
  enrichment = enrichment[order(enrichment$NES), ]
  idx = c()
  if(bottom>0) idx = 1:min(bottom, nrow(enrichment))
  if(top>0) idx = c(idx, max(1, nrow(enrichment)-top+1):nrow(enrichment))
  idx = unique(idx)
  enrichment = enrichment[idx, ]

  #normalize term description
  terms=as.character(enrichment$Description)
  terms=lapply(terms, function(x,k){
    x=as.character(x)
    if(nchar(x)>k){x=substr(x,start=1,stop=k)}
    return(x)}, charLength)
  enrichment$Description=do.call(rbind,terms)
  idx = duplicated(enrichment$Description)
  enrichment=enrichment[!idx,]
  enrichment$Name=factor(enrichment$Description, levels=enrichment$Description)

  #Prepare data for plotting.
  geneNames = strsplit(enrichment$geneName, "\\/")
  geneIds = strsplit(enrichment$geneID, "\\/")
  gg = data.frame(ID = rep(enrichment$ID, enrichment$Count), Term = rep(enrichment$Description, enrichment$Count),
             Gene = unlist(geneNames), geneIds = unlist(geneIds), stringsAsFactors = FALSE)
  if(keytype == "Symbol") gg$GeneScore = geneList[gg$Gene]
  if(keytype == "Entrez") gg$GeneScore = geneList[gg$geneIds]
  gg$Size = abs(gg$GeneScore)
  idx = is.na(gg$GeneScore) | (gg$GeneScore>gene_cutoff[1] & gg$GeneScore<gene_cutoff[2])
  gg = gg[!idx, ]
  gg = gg[order(gg$GeneScore), ]
  gg$Gene = factor(gg$Gene, levels = unique(gg$Gene))
  gg$Term = factor(gg$Term, levels = unique(gg$Term))
  # Plot the dot heatmap
  p1 = ggplot(data=gg, aes(x=Gene, y=Term, size=Size, color = GeneScore))
  p1 = p1 + geom_point()
  p1 = p1 + scale_color_gradient2(low = "#081087", high = "#c12603")
  p1 = p1 + theme(panel.grid.major=element_line(colour="gray90"),
                  panel.grid.minor=element_blank(),
                  panel.background=element_blank())
  p1 = p1 + labs(x=NULL, y=NULL, color = NULL)
  # p1 = p1 + theme(legend.position="top")
  p1 = p1 + scale_size_continuous(guide = FALSE)
  p1 = p1 + theme(legend.key = element_rect(fill = "transparent", colour = "transparent"))
  p1 = p1 + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                  plot.title = element_text(hjust = 0.5, size=18),
                  axis.text = element_text(colour="gray10"),
                  axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
  p1 = p1 + theme(axis.line = element_line(size=0.5, colour = "black"),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.border = element_blank(), panel.background = element_blank(),
                  legend.key = element_rect(fill = "transparent"))

  if(!is.null(filename)){
    ggsave(plot=p1, filename=filename, units = "in", width=width, height=height, ...)
  }
  return(p1)
}
