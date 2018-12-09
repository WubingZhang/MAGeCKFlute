#' View enriched terms
#'
#' Grid plot for enriched terms
#'
#' @docType methods
#' @name EnrichedView
#' @rdname EnrichedView
#' @aliases enrichview
#' @param enrichment A data frame of enrichment result, with columns of ID, Description, p.adjust and Count.
#' @param plotTitle Same as 'title' in 'plot'.
#' @param type "All" (Default), enriched geneset category for visualization.
#' @param color Color of nodes.
#' @param termNum Integer, specifying number of top enriched terms to show.
#' @param charLength Integer, specifying max length of enriched term name to show as coordinate lab.
#' @param filename Figure file name to create on disk. Default filename="NULL", which means
#' no output.
#' @param width As in ggsave.
#' @param height As in ggsave.
#' @param ... Other available parameters in ggsave.
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#' @author Feizhen Wu
#' @seealso \code{\link{KeggPathwayView}}
#' @seealso \code{\link{EnrichedGSEView}}
#' @examples
#' data(geneList, package = "DOSE")
#' enrichRes <- enrich.HGT(geneList[1:100])
#' EnrichedView(enrichment=enrichRes@result)
#' @import DOSE
#' @export

EnrichedView=function(enrichment, plotTitle = NULL, type = "All", color = "#3f90f7", termNum = 15,
                      charLength = 40, filename = NULL, width = 7, height = 4, ...){

  # Select subset of enrichedRes.
  if(type == "All") type = 'KEGG+BIOCARTA+REACTOME+GOBP+GOCC+GOMF+EHMN+PID+WikiPathways'
  type = unlist(strsplit(type, "\\+"))
  idx = toupper(gsub("_.*", "", enrichment$ID)) %in% toupper(type)
  enrichment = enrichment[idx, ]

  # No enriched pathways
  if(is.null(enrichment) || nrow(enrichment)==0){
    p1 = noEnrichPlot(plotTitle)
    if(!is.null(filename)){
      ggsave(plot=p1,filename=filename, units = "in", width=width, height=height, ...)
    }
    return(p1)
  }

  #The column of Description, ID, p.adjust, and Count are neccessary.
  enrichment$logP = -log10(enrichment$p.adjust)
  enrichment = enrichment[!is.na(enrichment$ID),]
  enrichment=enrichment[!duplicated(enrichment$Description),]
  enrichment = enrichment[order(enrichment$logP,decreasing = TRUE), ]
  if(nrow(enrichment)>=termNum){
    enrichment=enrichment[1:termNum,]
  }

  #normalize term description
  {
    terms=as.character(enrichment$Description)
    terms=lapply(terms,function(x,k){
      x=as.character(x)
      # x=paste(toupper(substr(x,1,1)),substr(x,2,nchar(x)),sep="")
      if(nchar(x)>k){
        x=substr(x,start=1,stop=k)
        # x=gsub("(\\w+)$", "...", x)
      }
      return(x)
    },charLength)
    enrichment$Description=do.call(rbind,terms)
  }
  idx = duplicated(enrichment$Description)
  enrichment=enrichment[!idx,]
  enrichment$Name=factor(enrichment$Description,levels=rev(enrichment$Description))

  p1 = ggplot(data=enrichment,aes(x=logP,y=Name,size=Count))
  p1 = p1 + geom_point(color=color)
  p1 = p1 + theme(panel.grid.major=element_line(colour="gray90"),
                 panel.grid.minor=element_blank(),
                 panel.background=element_blank())
  p1 = p1 + xlab(expression(-Log["10"](Adjust.pvalue)))+ylab("")
  p1 = p1 + labs(title=plotTitle)
  p1 = p1 + theme(legend.position="right")
  p1 = p1 + theme(legend.key = element_rect(fill = "transparent", colour = "transparent"))
  p1 = p1 + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p1 = p1 + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(),
                legend.key = element_rect(fill = "transparent"))

  if(!is.null(filename)){
    ggsave(plot=p1,filename=filename, units = "in", width=width, height=height, ...)
  }
  return(p1)
}


#' View enriched terms in GSEA
#'
#' Grid plot for enriched terms in GSEA
#'
#' @docType methods
#' @name EnrichedGSEView
#' @rdname EnrichedGSEView
#' @aliases enrichgseview
#' @param enrichment A data frame of enrichment result, with columns of ID, Description, p.adjust and NES
#' @param decreasing Logical. Should the NES sort order be increasing or decreasing?
#' @param plotTitle Same as 'title' in 'plot'.
#' @param type "All" (Default), enriched geneset category for visualization.
#' @param termNum Integer, specifying number of top enriched terms to show.
#' @param charLength Integer, specifying max length of enriched term name to show as coordinate lab.
#' @param filename Figure file name to create on disk. Default filename="NULL", which means
#' no output.
#' @param width As in ggsave.
#' @param height As in ggsave.
#' @param ... Other available parameters in ggsave.
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#' @author Wubing Zhang
#' @seealso \code{\link{EnrichedView}}
#' @examples
#' \dontrun{
#'     data(geneList, package = "DOSE")
#'     enrichRes = enrich.GSE(geneList, type = "KEGG", organism="hsa")
#'     EnrichedGSEView(as.data.frame(enrichRes), plotTitle = "GSEA Analysis")
#' }
#' @export
EnrichedGSEView = function(enrichment, decreasing = TRUE, plotTitle = NULL,
                         type = "All", termNum = 15, charLength = 40,
                         filename = NULL, width = 7, height = 4, ...){

  # Select subset of enrichedRes.
  if(type == "All") type = 'KEGG+BIOCARTA+REACTOME+GOBP+GOCC+GOMF+EHMN+PID+WikiPathways'
  type = unlist(strsplit(type, "\\+"))
  idx = toupper(gsub("_.*", "", enrichment$ID)) %in% toupper(type)
  enrichment = enrichment[idx, ]

  # No enriched pathways
  if(is.null(enrichment) || nrow(enrichment)==0){
    p1 = noEnrichPlot(plotTitle)
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
  enrichment = enrichment[order(enrichment$NES, decreasing = decreasing), ]
  if(nrow(enrichment) >= termNum){ enrichment = enrichment[1:termNum,] }

  # Normalize term description
  terms = as.character(enrichment$Description)
  terms = lapply(terms, function(x,k){
    x = as.character(x)
    x = paste(toupper(substr(x,1,1)),substr(x,2,nchar(x)), sep="")
    if(nchar(x)>k){ x = paste0(substr(x,start=1,stop=k), "...")}
    return(x)}, charLength)
  enrichment$Description = do.call(rbind, terms)
  idx = duplicated(enrichment$Description)
  enrichment=enrichment[!idx, ]
  enrichment$Name = factor(enrichment$Description, levels=rev(enrichment$Description))

  # Visualize the top enriched terms
  enrichment$col = "#e41a1c"
  enrichment$col[enrichment$NES<0] = "#377eb8"
  p1 <- ggplot(data=enrichment, aes(x=NES, y=Name, size=logP))
  p1 <- p1 + geom_point(color = enrichment$col)
  p1 <- p1 + theme(panel.grid.major=element_line(colour="gray90"),
                 panel.grid.minor=element_blank(),
                 panel.background=element_blank())
  p1 <- p1 + xlab("NES") + ylab(NULL)
  p1 <- p1 + labs(title=plotTitle)
  p1 = p1 + theme(legend.position="right")
  p1 = p1 + theme(legend.key = element_rect(fill = "transparent", colour = "transparent"))
  p1 = p1 + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p1 = p1 + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(),
                legend.key = element_rect(fill = "transparent"))

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
#' @param plotTitle Same as 'title' in 'plot'.
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#' @author Wubing Zhang
noEnrichPlot <- function(plotTitle){
  p1 = ggplot()
  p1 = p1 + geom_text(aes(x=0,y=0,label="No enriched terms"), size=6)
  p1 = p1 + labs(title=plotTitle)
  p1 = p1 + theme(plot.title = element_text(size=12))
  p1 = p1 + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                  plot.title = element_text(hjust = 0.5, size=18),
                  axis.text = element_text(colour="gray10"))
  p1 = p1 + theme(axis.line = element_line(size=0.5, colour = "black"),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.border = element_blank(), panel.background = element_blank())
  p1
}
