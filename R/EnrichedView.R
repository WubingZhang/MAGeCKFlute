#' View enriched terms
#'
#' Grid plot for enriched terms
#'
#' @docType methods
#' @name EnrichedView
#' @rdname EnrichedView
#' @param enrichment A data frame of enrichment result, with columns of ID, Description, p.adjust and NES.

#' @param rank_by "pvalue" or "NES", specifying the indices for ranking pathways.
#' @param mode 1 or 2.
#' @param subset A vector of pathway ids.
#' @param top An integer, specifying the number of upregulated terms to show.
#' @param bottom An integer, specifying the number of downregulated terms to show.

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
#' data(geneList, package = "DOSE")
#' \dontrun{
#'     enrichRes = enrich.GSE(geneList, organism="hsa")
#'     EnrichedView(enrichRes, top = 5, bottom = 5)
#' }
#' @export
EnrichedView = function(enrichment,
                        rank_by = "pvalue",
                        mode = 1,
                        subset = NULL,
                        top = 0, bottom = 0,
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
  flag = ifelse("NES"%in%colnames(enrichment), TRUE, FALSE)
  if(!flag) colnames(enrichment)[colnames(enrichment)=="Count"] = "NES"
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

  # Select pathways to show ##
  if(bottom>0){
    tmp = enrichment[enrichment$NES<0, ]
    subset = c(subset, tmp$ID[1:min(nrow(tmp), bottom)])
  }
  if(top>0){
    tmp = enrichment[enrichment$NES>0, ]
    subset = c(subset, tmp$ID[1:min(nrow(tmp), top)])
  }
  if(!is.null(subset)){
    idx = enrichment$ID %in% subset
    if(sum(idx)==0) return(noEnrichPlot("No eligible terms!!!"))
    enrichment = enrichment[idx, ]
  }
  # enrichment$ID = factor(enrichment$ID, levels=c(pid_neg, pid_pos))
  # enrichment = enrichment[order(enrichment$ID), ]
  enrichment$Description = factor(enrichment$Description,
                                  levels=unique(enrichment$Description))
  enrichment$ID = factor(enrichment$ID, levels=unique(enrichment$ID))

  ## Prepare data for plotting ##
  if(x=="NES"){
    enrichment$x = enrichment$NES
    enrichment$size = enrichment$logFDR
  }else if(x=="LogP"){
    enrichment$x = enrichment$logP
    enrichment$size = abs(enrichment$NES)
  }else{
    enrichment$x = enrichment$logFDR
    enrichment$size = abs(enrichment$NES)
  }

  # Visualize the top enriched terms
  enrichment$col = "UG"
  enrichment$col[enrichment$NES<0] = "DG"

  if(mode==1){
    ## Plot the figure ##
    p1 = ggplot(data=enrichment, aes_string(x="x", y="Description", size="size"))
    p1 = p1 + geom_point(aes(color = col))
    p1 = p1 + scale_color_manual(values = c("DG"="#377eb8", "UG"="#e41a1c"))
    p1 = p1 + theme(panel.grid.major=element_line(colour="gray90"),
                    panel.grid.minor=element_blank(),
                    panel.background=element_blank())
    p1 = p1 + theme(legend.position="right")
    p1 = p1 + theme(legend.key = element_rect(fill = "transparent", colour = "transparent"))
    p1 = p1 + theme_bw(base_size = 14)
    p1 = p1 + theme(plot.title = element_text(hjust = 0.5))
    p1 = p1 + guides(color = FALSE)
  }else if(mode == 2){
    idx = (max(enrichment$x)-enrichment$x) > enrichment$x
    enrichment$hjust = 1.1
    enrichment$hjust[idx] = -0.1
    p1 = ggplot(enrichment, aes_string("x", "ID", label = "Description"))
    p1 = p1 + geom_point(aes_string(color = "col", size = "size"))
    p1 = p1 + xlim(0, NA)
    p1 = p1 + geom_text(aes_string(hjust = "hjust"))
    p1 = p1 + theme_bw(base_size = 14) + theme(plot.title = element_text(hjust = 0.5))
  }
  if(x=="NES"){
    p1 = p1 + labs(x = "Enrichment score", y = NULL, color = NULL,
                   size = expression(-log[10]*FDR))
    if(!flag) p1 = p1 + labs(x = "Count")
  }else{
    p1 = p1 + labs(x = expression(-log[10]*FDR), y = NULL, color = NULL, size = "NES")
    if(!flag) p1 = p1 + labs(size = "Count")
  }

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
  p1 = p1 + theme_bw(base_size = 14)
  p1 = p1 + theme(plot.title = element_text(hjust = 0.5))
  p1
}
