#' Violin plot
#'
#' Plots the violin of beta scores in Control and Treatment samples.
#'
#' @docType methods
#' @name ViolinView
#' @rdname ViolinView
#' @aliases violinview
#'
#' @param beta data frame, which has columns of 'Gene', \code{samples}.
#' @param samples character, specifying the name of samples to be compared
#' @param main as in 'plot'.
#' @param ylab as in 'plot'.
#' @param filename figure file name to create on disk. Default filename="NULL", which means
#' don't save the figure on disk.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#' @note See the vignette for an example of density plot for beta score deviation.
#' Note that the source code of \code{ViolinView} is very simple.
#' The source can be found by typing \code{MAGeCKFlute:::ViolinView}
#' or \code{getMethod("ViolinView")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/ViolinView.R}
#' Users should find it easy to customize this function.
#'
#' @seealso \code{\link{DensityDiffView}}   \code{\link{DensityView}}
#' @seealso \code{\link{MAView}}   \code{\link{SquareView}}
#' @seealso \code{\link{CellCycleView}}  \code{\link{EnrichedView}}
#' @seealso \code{\link{EnrichedGSEView}}  \code{\link{KeggPathwayView}}
#' @seealso \code{\link{RankView}}    \code{\link{ScatterView}}
#'
#' @examples
#' data(MLE_Data)
#' # Read beta score from gene summary table in MAGeCK MLE results
#' dd = ReadBeta(MLE_Data, organism="hsa")
#' dd = dd[,-2]
#' ViolinView(dd)
#'
#'
#' @importFrom reshape melt
#' @importFrom ggsci scale_color_npg
#'
#' @export
#'

#===Distribution of beta scores======================================
ViolinView <- function(beta, samples=NULL, main=NULL,ylab="Beta Score",filename=NULL){
  requireNamespace("reshape", quietly=TRUE) || stop("need reshape package")
  requireNamespace("ggsci", quietly=TRUE) || stop("need ggsci package")

  dd1=beta
  loginfo(paste("Violin plot for", main, ylab, "..."))
  if(!is.null(samples) && length(samples)>1){ dd1 = dd1[,c("Gene", samples)]}

  dd1 = melt(dd1,id="Gene")
  #======
  p=ggplot(data=dd1,aes(x=variable,y=value,color=variable))
  p=p+geom_violin()+geom_boxplot(width=.1, outlier.colour=NA)
  #p=p+ylim(-1.5,1)
  p=p+scale_color_npg()+theme_bw(14)+theme(plot.title = element_text(hjust = 0.5,size=12))
  p=p+theme(legend.position = "none")
  p=p+theme(axis.text.x = element_text(angle = 30, hjust = 1, size=8))
  p=p+labs(x="Sample",y=ylab,title=main)

  if(!is.null(filename)){
    ggsave(plot=p,filename=filename,units = "in",width=300/100,height =270/100 )
  }
  return(p)
}
