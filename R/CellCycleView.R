#' Plot of linear fitting lines for beta scores in Control and Treatment sample
#'
#' Estimate cell cycle time in different samples by linear fitting of beta scores, and plot fitting lines,
#' in which x-axis is control beta score and y-axis is beta score of all samples.
#'
#' @docType methods
#' @name CellCycleView
#' @rdname CellCycleView
#' @aliases CellCycle,MAGeCKFlute-method
#'
#' @param beta data frame, which has columns of 'Gene', \code{ctrlname} and other samples.
#' @param ctrlname character vector, specifying the name of control sample.
#' @param main as in 'plot'.
#' @param filename figure file name to create on disk. Default filename="NULL", which means
#' no output.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#' @note See the vignette for an example of cell cycle estimate and view.
#' Note that the source code of \code{CellCycleFit} is very simple.
#' The source can be found by typing \code{MAGeCKFlute:::CellCycleView}
#' or \code{getMethod("CellCycleView")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/CellCycleView.R}
#' Users should find it easy to customize this function.
#'
#' @seealso \code{\link{DensityDiffView}}   \code{\link{DensityView}}
#' @seealso \code{\link{ViolinView}}   \code{\link{SquareView}}
#' @seealso \code{\link{KeggPathwayView}}  \code{\link{EnrichedView}}
#' @seealso \code{\link{EnrichedGSEView}}  \code{\link{MAView}}
#' @seealso \code{\link{RankView}}    \code{\link{ScatterView}}
#'
#' @examples
#' data(MLE_Data)
#' # Read beta score from gene summary table in MAGeCK MLE results
#' dd = ReadBeta(MLE_Data, organism="hsa")[,-2]
#' CellCycleView(dd, ctrlname = c("D7_R1", "D7_R2"))
#'
#' @importFrom reshape melt
#'
#' @export

#===Distribution of beta scores======================================
CellCycleView <- function(beta, ctrlname="Control", main=NULL, filename=NULL){
  loginfo(paste("Cell cycle fitting of treatments compaired to",
                paste(ctrlname, collapse = "&")))
  # dd2 = beta[,c("Gene", ctrlname, treatname)]
  requireNamespace("reshape", quietly=TRUE) || stop("need reshape package")
  dd2 = beta
  idx = colnames(beta) %in% ctrlname
  dd2 = melt(dd2[,!idx], id="Gene")
  dd2$x = rep(rowMeans(beta[,ctrlname,drop=FALSE]), nrow(dd2)/nrow(beta))

  p=ggplot(dd2,aes(x,value,color=variable,group=variable))
  p=p+geom_point(alpha=4/10,size=0.8)
  p=p+geom_smooth(method='lm',se=FALSE)
  p=p+geom_abline(slope = 1, intercept = 0, color="gray50", linetype=2)
  p=p+theme_bw(12)+theme(plot.title = element_text(hjust = 0.5,size=12))
  p=p+labs(x="Control", y="Treatment", title=main, color=NULL)
  p=p+theme(legend.justification = c(0, 1), legend.position = c(0.01, 0.99))

  if(!is.null(filename)){
    ggsave(plot=p,filename=filename,units = "in",width=300/100,height =270/100 )
  }
  return(p)
}


