#' Plot rank of gene points
#'
#' Rank all genes according to beta score deviation, and label top and bottom meaningful genes.
#' Some other interested genes can be labeled too.
#'
#' @docType methods
#' @name RankView
#' @rdname RankView
#' @aliases rankview
#'
#' @param beta data frame containing columns of "Gene" and "diff".
#' @param genelist character vector, specifying labeled genes besides top and bottom labeled genes.
#' @param top integer, specifying top number of genes to be labeled.
#' @param bottom integer, specifying bottom number of genes to be labeled.
#' @param cutoff numeric, cutoff of \code{diff}
#' @param main as in 'plot'.
#' @param filename figure file name to create on disk. Default filename="NULL", which means
#' no output.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#' @note See the vignette for an example of CellCycleView.
#' Note that the source code of \code{RankView} is very simple.
#' The source can be found by typing \code{MAGeCKFlute:::RankView}
#' or \code{getMethod("RankView")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/CellCycleView.R}
#' Users should find it easy to customize this function.
#'
#' @seealso \code{\link{DensityDiffView}}   \code{\link{DensityView}}
#' @seealso \code{\link{ViolinView}}   \code{\link{SquareView}}
#' @seealso \code{\link{CellCycleView}}  \code{\link{EnrichedView}}
#' @seealso \code{\link{EnrichedGSEView}}  \code{\link{KeggPathwayView}}
#' @seealso \code{\link{MAView}}    \code{\link{ScatterView}}
#'
#' @examples
#' data(MLE_Data)
#' # Read beta score from gene summary table in MAGeCK MLE results
#' dd = ReadBeta(MLE_Data, ctrlName = "D7_R1", treatName = "PLX7_R1", organism="hsa")
#' dd$diff = dd$Treatment - dd$Control
#' RankView(dd[,c("Gene", "diff")])
#'
#' @importFrom reshape melt
#' @importFrom ggrepel geom_label_repel
#'
#' @export
#'

RankView <- function(beta, genelist=c(), top=10, bottom=10,
                     cutoff=c(sd(beta$diff), -sd(beta$diff)), main=NULL,filename=NULL){
  requireNamespace("ggrepel", quietly=TRUE) || stop("need ggrepel package")
  requireNamespace("reshape", quietly=TRUE) || stop("need reshape package")
  loginfo(paste("Rank genes and plot..."))
  diff = as.numeric(beta$diff)
  data = list()
  data$diff = diff
  data$Rank = rank(diff)
  data$Gene = beta$Gene
  data$group = "no"
  data = as.data.frame(data, stringsAsFactors=FALSE)
  data$group[data$diff>cutoff[1]] = "up"
  data$group[data$diff<cutoff[2]] = "down"


  idx=(data$Rank<=bottom) | (data$Rank>(max(data$Rank)-top)) | (data$Gene %in% genelist)
  mycolour=c("no"="darkgray",  "up"="#e41a1c","down"="#377eb8")
  p=ggplot()
  p=p+geom_point(aes(x=diff,y=Rank,color=factor(data$group)),data=data,size=0.5)
  p=p+theme_bw(14)+theme(plot.title = element_text(hjust = 0.5,size=12))
  p=p+scale_color_manual(values=mycolour)
  p=p+geom_jitter(width = 0.1, height = 0.1)
  p=p+theme(panel.background = element_rect(fill='white', colour='black'))
  p=p+geom_vline(xintercept = 0,linetype = "dotted")+geom_vline(xintercept = cutoff,linetype = "dotted")
  p=p+geom_label_repel(aes(x=diff, y=Rank,fill=group,label = Gene),data=data[idx,],
                       fontface = 'bold', color = 'white',size = 2.5,force=3,
                       box.padding = unit(0.4, "lines"),segment.color = 'grey50',
                       point.padding = unit(0.3, "lines"))
  p=p+scale_fill_manual(values=mycolour)
  p=p+labs(x="Treatment-Control Beta Score",y="Rank",title=main)
  p=p+theme(legend.position="none")#+ylim(-1000,7000)
  if(!is.null(filename)){
    ggsave(plot=p,filename=filename,units = "in",width=5,height=4 )
  }
  return(p)
}
