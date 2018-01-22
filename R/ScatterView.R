#' Scatter plot
#'
#' Scatter plot of all genes, in which x-axis is mean beta score in Control samples, y-axis
#' is mean beta scores in Treatment samples.
#'
#' @docType methods
#' @name ScatterView
#' @rdname ScatterView
#' @aliases scatterview
#'
#' @param beta Data frame, which has columns of 'Gene', \code{ctrlname} and \code{treatname}.
#' @param ctrlname A character, specifying the name of control sample.
#' @param treatname A character, specifying the name of treatment sample.
#' @param scale_cutoff Boolean or numeric, whether scale cutoff to whole genome level,
#' or how many standard deviation will be used as cutoff.
#' @param main As in 'plot'.
#' @param filename Figure file name to create on disk. Default filename="NULL", which means
#' don't save the figure on disk.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#' @note See the vignette for an example of ScatterView.
#' Note that the source code of \code{ScatterView} is very simple.
#' The source can be found by typing \code{MAGeCKFlute:::ScatterView}
#' or \code{getMethod("ScatterView")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/ScatterView.R}
#' Users should find it easy to customize this function.
#'
#' @seealso \code{\link{DensityDiffView}}   \code{\link{DensityView}}
#' @seealso \code{\link{ViolinView}}   \code{\link{SquareView}}
#' @seealso \code{\link{CellCycleView}}  \code{\link{EnrichedView}}
#' @seealso \code{\link{EnrichedGSEView}}  \code{\link{KeggPathwayView}}
#' @seealso \code{\link{RankView}}    \code{\link{MAView}}
#'
#' @examples
#' data(MLE_Data)
#' # Read beta score from gene summary table in MAGeCK MLE results
#' dd = ReadBeta(MLE_Data, organism="hsa")
#' ScatterView(dd, ctrlname = "D7_R1", treatname = "PLX7_R1")
#'
#'
#' @export
#'

ScatterView <- function(beta, ctrlname="Control",treatname="Treatment", scale_cutoff=1,
                     main=NULL,filename=NULL){

  beta$Control=rowMeans(beta[,ctrlname,drop= FALSE])
  beta$Treatment=rowMeans(beta[,treatname,drop= FALSE])
  intercept=CutoffCalling(beta$Treatment-beta$Control, scale=scale_cutoff)
  beta$diff = beta$Treatment - beta$Control
  beta$group="no"
  beta$group[beta$diff>intercept]="up"
  beta$group[beta$diff<(-intercept)]="down"

  data=beta
  loginfo(paste("Scatter plot of", main, "Treat-Ctrl beta scores ..."))
  mycolour=c("no"="aliceblue",  "up"="#e41a1c","down"="#377eb8")
  xmin=min(data$Control)
  xmax=max(data$Control)
  ymin=min(data$Treatment)
  ymax=max(data$Treatment)
  #=========
  p=ggplot(data,aes(x=Control,y=Treatment,colour=group,fill=group))
  p=p+geom_point(position = "identity",shape=".",alpha=1/100,size = 0.01,show.legend = FALSE)
  p=p+scale_color_manual(values=mycolour)
  p=p+geom_jitter(position = "jitter",show.legend = FALSE)
  p=p+theme_bw(14)+theme(plot.title = element_text(hjust = 0.5,size=12))
  p=p+geom_abline(intercept = -intercept)
  p=p+geom_abline(intercept = +intercept)
  p=p+labs(x="Control beta.score",y="Treatment beta.score",title=main)
  p=p+annotate("text",color="#e41a1c",x=xmin, y=ymax,hjust = 0,
               label=paste("GroupA: ",as.character(dim(data[data$group=="up",])[1]),sep=""))
  p=p+annotate("text",color="#377eb8",x=xmax, y=ymin,hjust = 1,
               label=paste("GroupB: ",as.character(dim(data[data$group=="down",])[1]),sep=""))
  #============
  if(!is.null(filename)){
    write.table(beta, file.path(dirname(filename), paste0("GroupAB_", main, ".txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    ggsave(plot=p,filename=filename,units = "in",width=5,height =4 )
  }
  return(p)
}

