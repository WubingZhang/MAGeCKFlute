#' View enriched terms
#'
#' Grid plot for enriched terms
#'
#' @docType methods
#' @name EnrichedView
#' @rdname EnrichedView
#' @param enrichment A data frame of enrichment result, with columns of ID, Description, p.adjust and NES.

#' @param rank_by "pvalue" or "NES", specifying the indices for ranking pathways.
#' @param top An integer, specifying the number of top enriched terms to show.
#' @param bottom An integer, specifying the number of bottom enriched terms to show.

#' @param x Character, "NES", "LogP", or "LogFDR", indicating the variable on the x-axis.
#' @param charLength Integer, specifying max length of enriched term name to show as coordinate lab.

#' @param filename Figure file name to create on disk. Default filename="NULL".
#' @param width As in ggsave.
#' @param height As in ggsave.
#' @param ... Other available parameters in ggsave.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#' @author Wubing Zhang
#' @seealso \code{\link{EnrichedView}}
#' @examples
#' \dontrun{
#'     data(geneList, package = "DOSE")
#'     enrichRes = enrich.GSE(geneList, organism="hsa")
#'     EnrichedView(slot(enrichRes, "result"))
#' }
#' @export
EnrichedView = function(enrichment,
                        rank_by = "pvalue",
                        top = 5, bottom = 0,
                        x = "LogFDR",
                        charLength = 40,
                        filename = NULL,
                        width = 7, height = 4, ...){

  if(is(enrichment, "enrichResult")) enrichment = slot(enrichment, "result")
  if(is(enrichment, "gseaResult")) enrichment = slot(enrichment, "result")
  # No enriched pathways
  if(is.null(enrichment) || nrow(enrichment)==0){
    p1 = noEnrichPlot("No enriched terms")
    if(!is.null(filename)){
      ggsave(plot=p1,filename=filename, units = "in", width=width, height=height, ...)
    }
    return(p1)
  }

  ## Rank enriched pathways ##
  enrichment$logP = round(-log10(enrichment$pvalue), 1)
  enrichment$logFDR = round(-log10(enrichment$p.adjust), 1)
  enrichment = enrichment[!is.na(enrichment$ID), ]
  if(tolower(rank_by) == "pvalue"){
    enrichment = enrichment[order(enrichment$pvalue, -abs(enrichment$NES)), ]
  }else if(tolower(rank_by) == "nes"){
    enrichment = enrichment[order(-abs(enrichment$NES), enrichment$pvalue), ]
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
  idx = enrichment$ID %in% c(pid_neg, pid_pos)
  if(sum(idx)==0) return(noEnrichPlot("No eligible terms!!!"))
  enrichment = enrichment[idx, ]
  enrichment$ID = factor(enrichment$ID, levels=c(pid_neg, pid_pos))
  enrichment = enrichment[order(enrichment$ID), ]
  enrichment$Description = factor(enrichment$Description, levels=unique(enrichment$Description))

  ## Prepare data for plotting ##
  if(x=="NES"){
    enrichment$x = enrichment$NES
    enrichment$size = enrichment$logFDR
  }else if(x=="LogP"){
    enrichment$x = enrichment$logP
    enrichment$size = enrichment$NES
  }else{
    enrichment$x = enrichment$logFDR
    enrichment$size = enrichment$NES
  }

  # Visualize the top enriched terms
  enrichment$col = "Up"
  enrichment$col[enrichment$NES<0] = "Down"

  ## Plot the figure ##
  p1 = ggplot(data=enrichment, aes(x=x, y=Description, size=size))
  p1 = p1 + geom_point(aes(color = col))
  p1 = p1 + scale_color_manual(values = c("Down"="#377eb8", "Up"="#e41a1c"))
  p1 = p1 + theme(panel.grid.major=element_line(colour="gray90"),
                   panel.grid.minor=element_blank(),
                   panel.background=element_blank())
  p1 = p1 + theme(legend.position="right")
  p1 = p1 + theme(legend.key = element_rect(fill = "transparent", colour = "transparent"))
  p1 = p1 + theme(text = element_text(colour="black",size = 11, family = "Helvetica"),
                  plot.title = element_text(hjust = 0.5, size=14),
                  axis.text = element_text(colour="gray10"))
  p1 = p1 + theme(axis.line = element_line(size=0.5, colour = "black"),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.border = element_blank(), panel.background = element_blank(),
                  legend.key = element_rect(fill = "transparent"))
  p1 = p1 + guides(color = FALSE)
  if(x=="NES")
    p1 = p1 + labs(x = "Enrichment score", y = NULL, col = "Group", size = "LogFDR")
  else
    p1 = p1 + labs(x = x, y = NULL, col = "Group", size = "NES")

  if(!is.null(filename)){
    ggsave(plot=p1, filename=filename, units = "in", width=width, height=height, ...)
  }
  return(p1)
}

#' Blank figure
#'
#' @docType methods
#' @name noEnrichPlot
#' @rdname noEnrichPlot
#' @param main The title of figure.
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#' @author Wubing Zhang
noEnrichPlot = function(main = "No enriched terms"){
  p1 = ggplot()
  p1 = p1 + geom_text(aes(x=0, y=0, label="No enriched terms"), size=6)
  p1 = p1 + labs(title=main)
  p1 = p1 + theme(plot.title = element_text(size=12))
  p1 = p1 + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                  plot.title = element_text(hjust = 0.5, size=18),
                  axis.text = element_text(colour="gray10"))
  p1 = p1 + theme(axis.line = element_line(size=0.5, colour = "black"),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.border = element_blank(), panel.background = element_blank())
  p1
}
