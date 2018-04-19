#' Estimate cell cycle time for all samples compared to control sample and view.
#'
#' Estimate cell cycle time in different samples by linear fitting of beta scores, and plot fitting lines,
#' in which x-axis is control beta score and y-axis is beta score of all samples.
#'
#' @docType methods
#' @name CellCycleView
#' @rdname CellCycleView
#' @aliases CellCycle,MAGeCKFlute-method
#'
#' @param beta Data frame, which has columns of 'Gene', \code{ctrlname} and other samples.
#' @param ctrlname Character vector, specifying the name of control sample.
#' @param main As in 'plot'.
#' @param filename Figure file name to create on disk. Default filename="NULL", which means
#' no output.
#' @param width As in ggsave.
#' @param height As in ggsave.
#' @param ... Other available parameters in ggsave.
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
#'
#' @examples
#' data(MLE_Data)
#' # Read beta score from gene summary table in MAGeCK MLE results
#' dd = ReadBeta(MLE_Data, organism="hsa")
#' CellCycleView(dd[,-2], ctrlname = c("D7_R1", "D7_R2"))
#'
#' @importFrom reshape melt
#'
#' @export

#===Distribution of beta scores======================================
CellCycleView <- function(beta, ctrlname="Control", main=NULL, filename=NULL, width=5, height = 4, ...){
  loginfo(paste("Cell cycle fitting of treatments compaired to",
                paste(ctrlname, collapse = "&")))
  # dd2 = beta[,c("Gene", ctrlname, treatname)]
  requireNamespace("reshape", quietly=TRUE) || stop("need reshape package")
  dd2 = beta
  idx = colnames(beta) %in% ctrlname
  dd2 = melt(dd2[,!idx], id="Gene")
  dd2$x = rep(rowMeans(beta[,ctrlname,drop=FALSE]), nrow(dd2)/nrow(beta))

  p = ggplot(dd2,aes(x,value,color=variable,group=variable))
  p = p + geom_point(alpha=4/10,size=0.5)
  p = p + geom_smooth(method='lm',se=FALSE, size=0.5)
  p = p + geom_abline(slope = 1, intercept = 0, color="gray50", linetype=2)
  p = p + labs(x="Control", y="Treatment", title=main, color=NULL)
  p = p + theme(legend.key = element_rect(fill = "transparent"))
  p = p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank())
  p = p + theme(legend.justification = c(0, 1), legend.position = c(0.01, 0.99))

  if(!is.null(filename)){
    ggsave(plot=p,filename=filename, units = "in", dpi=600, width=width, height = height, ...)
  }
  return(p)
}


