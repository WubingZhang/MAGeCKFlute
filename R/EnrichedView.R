#' View enriched terms
#'
#' Grid plot for enriched terms
#'
#' @docType methods
#' @name EnrichedView
#' @rdname EnrichedView
#' @aliases enrichview
#'
#' @param enrichment A data frame of enrichment result, with columns of ID, Description, p.adjust and Count.
#' @param plotTitle Same as 'title' in 'plot'.
#' @param color Color of nodes.
#' @param termNum Integer, specifying number of top enriched terms to show.
#' @param charLength Integer, specifying max length of enriched term name to show as coordinate lab.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Feizhen Wu
#'
#' @note  See the vignette for an example of EnrichedView
#' The source can be found by typing \code{MAGeCKFlute:::EnrichedView}
#' or \code{getMethod("EnrichedView")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/EnrichedView.R}
#' Users should find it easy to customize this function.
#'
#' @seealso \code{\link{KeggPathwayView}}
#' @seealso \code{\link{EnrichedGSEView}}
#'
#' @examples
#' data(MLE_Data)
#' universe = TransGeneID(MLE_Data$Gene, "SYMBOL", "ENTREZID", organism = "hsa")
#' genes = universe[1:50]
#' enrichRes <- enrich.HGT(genes, universe)
#' EnrichedView(enrichment=enrichRes@result, plotTitle = "Hypergemetric test")
#'
#'
#' @export


EnrichedView=function(enrichment,plotTitle=NULL,color="#3f90f7",termNum=20,charLength=40){

  if(is.null(enrichment) || nrow(enrichment)==0){
    p1=ggplot()
    p1=p1+geom_text(aes(x=0,y=0,label="No enriched terms"),size=6)
    p1=p1+labs(title=plotTitle)
    p1=p1+theme_bw(14)+theme(plot.title = element_text(hjust = 0.5,size=12))
    return(p1)
  }

  if(nrow(enrichment)>=termNum){
    enrichment=enrichment[1:termNum,]
  }

  #The column of Description, ID, p.adjust, and Count are neccessary.
  enrichment$logP = -log10(enrichment$p.adjust)
  enrichment = enrichment[!is.na(enrichment$ID),]
  enrichment=enrichment[!duplicated(enrichment$Description),]
  enrichment = enrichment[order(enrichment$logP,decreasing = TRUE), ]

  #normalize term description
  {
    terms=as.character(enrichment$Description)
    terms=lapply(terms,function(x,k){
      x=as.character(x)
      x=paste(toupper(substr(x,1,1)),substr(x,2,nchar(x)),sep="")
      if(nchar(x)>k){
        x=substr(x,start=1,stop=k)
        x=gsub("(\\w+)$", "...", x)
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
  p1 <- p1 + theme(panel.grid.major=element_line(colour="gray90"),
                 panel.grid.minor=element_blank(),
                 panel.background=element_blank())
  p1 <- p1 + xlab("-log10(Adjust.pvalue)")+ylab("")
  p1 <- p1 + labs(title=plotTitle)
  p1 <- p1 + theme_bw(8)
  # p1 <- p1 + theme(axis.text.x=element_text(size=8, face="plain", colour='black'))
  # p1 <- p1 + theme(axis.text.y=element_text(size=8, face="plain", colour='black'))
  p1 = p1 + theme(legend.position="bottom")
  p1 = p1 + theme_bw(14) + theme(plot.title = element_text(hjust = 0.5,size=12))

  return(p1)
}




##=====================================================================
#' View enriched terms in GSEA
#'
#' Grid plot for enriched terms in GSEA
#'
#' @docType methods
#' @name EnrichedGSEView
#' @rdname EnrichedGSEView
#' @aliases enrichgseview
#'
#' @param enrichment A data frame of enrichment result, with columns of ID, Description, p.adjust and NES
#' @param plotTitle Same as 'title' in 'plot'.
#' @param color Color of nodes
#' @param termNum Integer, specifying number of top enriched terms to show
#' @param charLength Integer, specifying max length of enriched term name to show as coordinate lab
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#' @note  See the vignette for an example of EnrichedGSEView
#' The source can be found by typing \code{MAGeCKFlute:::EnrichedGSEView}
#' or \code{getMethod("EnrichedGSEView")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/EnrichedGSEView.R}
#' Users should find it easy to customize this function.
#'
#' @seealso \code{\link{EnrichedView}}
#'
#' @examples
#' data(MLE_Data)
#' universe = TransGeneID(MLE_Data$Gene, "SYMBOL", "ENTREZID", organism = "hsa")
#' geneList = MLE_Data$D7_R1.beta
#' names(geneList) = universe
#' geneList = geneList[!is.na(universe)]
#' enrichRes = enrich.GSE(geneList, type = "KEGG", organism="hsa")
#' EnrichedGSEView(enrichRes@result, plotTitle = "GSEA Analysis")
#'
#'
#'
#' @export

##===================
EnrichedGSEView=function(enrichment,plotTitle=NULL, color="#3f90f7",termNum=20,charLength=50){

  if(nrow(enrichment)==0){
    p1=ggplot()
    p1=p1+geom_text(aes(x=0,y=0,label="No enriched terms"),size=6)
    p1=p1+labs(title=plotTitle)
    p1=p1+theme(plot.title = element_text(size=12))
    p1=p1+theme_bw(14)+theme(plot.title = element_text(hjust = 0.5,size=12))
    return(p1)
  }

  if(nrow(enrichment)>=termNum){
    enrichment=enrichment[1:termNum,]
  }

  #The column of Description, ID, p.adjust, and Count are neccessary.
  enrichment$logP = -log10(enrichment$p.adjust)
  enrichment = enrichment[!is.na(enrichment$ID),]
  enrichment=enrichment[!duplicated(enrichment$Description),]
  enrichment = enrichment[order(enrichment$logP,decreasing = TRUE), ]

  #normalize term description
  {
    terms=as.character(enrichment$Description)
    terms=lapply(terms,function(x,k){
      x=as.character(x)
      x=paste(toupper(substr(x,1,1)),substr(x,2,nchar(x)),sep="")
      if(nchar(x)>k){
        x=substr(x,start=1,stop=k)
        x=gsub("(\\w+)$", "...", x)
      }
      return(x)
    },charLength)
    enrichment$Description=do.call(rbind,terms)
  }
  idx = duplicated(enrichment$Description)
  enrichment=enrichment[!idx,]
  enrichment$Name=factor(enrichment$Description,levels=rev(enrichment$Description))

  p1=ggplot(data=enrichment,aes(x=logP,y=Name,size=NES))
  p1=p1+geom_point(color=color)
  p1 <- p1+theme(panel.grid.major=element_line(colour="gray90"),
                 panel.grid.minor=element_blank(),
                 panel.background=element_blank())
  p1 <- p1 + xlab("-log10(Adjust.pvalue)")+ylab("")
  p1 <- p1 + labs(title=plotTitle)
  p1 <- p1 + theme_bw(8)
  p1 = p1 + theme(legend.position="bottom")
  p1 = p1 + theme_bw(14) + theme(plot.title = element_text(hjust = 0.5,size=12))
  #p1
  return(p1)
}


