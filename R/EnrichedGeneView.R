#' Visualize enriched pathways and genes in those pathways
#'
#' @docType methods
#' @name EnrichedGeneView
#' @rdname EnrichedGeneView
#'
#' @param enrichment A data frame of enrichment result or an \code{enrichResult} object.
#' @param geneList A numeric geneList used in enrichment anlaysis.
#' @param rank_by "p.adjust" or "NES", specifying the indices for ranking pathways.
#'
#' @param top An integer, specifying the number of positively enriched terms to show.
#' @param bottom An integer, specifying the number of negatively enriched terms to show.
#' @param custom_pid A character vector (pathway IDs), customizing pathways to show.
#'
#' @param keytype "Entrez" or "Symbol".
#' @param gene_cutoff A two-length numeric vector, specifying cutoff for genes to show.
#' @param custom_gene A character vector (gene names), customizing genes to show.
#'
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
#' EnrichedGeneView(enrichment=slot(enrichRes, "result"), geneList, keytype = "Entrez")
#' @export

EnrichedGeneView=function(enrichment, geneList,
                          rank_by = "p.adjust",
                          top = 5, bottom = 5,
                          custom_pid = NULL,
                          keytype = "Symbol",
                          gene_cutoff = c(-log2(1.5), log2(1.5)),
                          custom_gene = NULL,
                          charLength = 40,
                          filename = NULL,
                          width = 7, height = 5, ...){

  if(is(enrichment, "enrichResult")) enrichment = enrichment@result
  if(is(enrichment, "gseaResult")) enrichment = enrichment@result

  ## No enriched pathways ##
  if(is.null(enrichment) || nrow(enrichment)==0){
    p1 = noEnrichPlot("No enriched terms")
    if(!is.null(filename)){
      ggsave(plot=p1,filename=filename, units = "in", width=width, height=height, ...)
    }
    return(p1)
  }

  ## Rank enriched pathways ##
  enrichment$logP = round(-log10(enrichment$p.adjust), 1)
  enrichment = enrichment[!is.na(enrichment$ID), ]
  if(tolower(rank_by) == "p.adjust"){
    enrichment = enrichment[order(enrichment$p.adjust), ]
  }else if(tolower(rank_by) == "nes"){
    enrichment = enrichment[order(abs(enrichment$NES), decreasing = TRUE), ]
  }

  ## Normalize term description ##
  terms = as.character(enrichment$Description)
  terms = lapply(terms, function(x,k){
    x = as.character(x)
    if(nchar(x)>k){x=substr(x,start=1,stop=k)}
    return(x)}, charLength)
  enrichment$Description = do.call(rbind, terms)
  enrichment = enrichment[!duplicated(enrichment$Description),]

  ## Select pathways to show ##
  pid_neg <- pid_pos <- NULL
  if(bottom>0){
    tmp = enrichment[enrichment$NES<0, ]
    pid_neg = tmp$ID[1:min(nrow(tmp), bottom)]
  }
  if(top>0){
    tmp = enrichment[enrichment$NES>0, ]
    pid_pos = tmp$ID[1:min(nrow(tmp), top)]
  }
  idx = enrichment$ID %in% c(custom_pid, pid_neg, pid_pos)
  if(sum(idx)==0) return(noEnrichPlot("No eligible terms!!!"))

  ## Prepare data for plotting ##
  enrichment = enrichment[idx, ]
  enrichment$Description = factor(enrichment$Description,
                                  levels=enrichment$Description)
  geneNames = strsplit(enrichment$geneName, "\\/")
  geneIds = strsplit(enrichment$geneID, "\\/")
  gg = data.frame(ID = rep(enrichment$ID, enrichment$Count),
                  Term = rep(enrichment$Description, enrichment$Count),
                  Size = rep(enrichment$logP, enrichment$Count),
                  Gene = unlist(geneNames), geneIds = unlist(geneIds),
                  stringsAsFactors = FALSE)

  ## Select genes to show ##
  names(geneList) = names(geneList)
  geneList = geneList[geneList<gene_cutoff[1] | geneList>gene_cutoff[2] |
                        names(geneList) %in% custom_gene]
  if(keytype == "Symbol") gg$GeneScore = geneList[gg$Gene]
  if(keytype == "Entrez") gg$GeneScore = geneList[gg$geneIds]

  ## Rank pathways and genes ##
  gg = gg[!is.na(gg$GeneScore), ]
  gg$Term = factor(gg$Term, levels = unique(gg$Term))
  gg = gg[order(gg$GeneScore), ]
  gg$Gene = factor(gg$Gene, levels = unique(gg$Gene))
  # Plot the dot heatmap
  p1 = ggplot(data=gg, aes(x=Gene, y=Term, size=Size, color = GeneScore))
  p1 = p1 + geom_point()
  p1 = p1 + scale_color_gradient2(low = "#081087", high = "#c12603")
  p1 = p1 + theme(panel.grid.major=element_line(colour="gray90"),
                  panel.grid.minor=element_blank(),
                  panel.background=element_blank())
  p1 = p1 + labs(x=NULL, y=NULL, color = "Gene score", size = "LogP")
  # p1 = p1 + theme(legend.position="top")
  # p1 = p1 + scale_size_continuous(guide = FALSE)
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
