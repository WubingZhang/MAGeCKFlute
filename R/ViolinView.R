#' Violin plot
#'
#' Plots the violin of beta scores in Control and Treatment samples.
#'
#' @docType methods
#' @name ViolinView
#' @rdname ViolinView
#' @aliases violinview
#'
#' @param beta Data frame, , including all \code{samples} as columns.
#' @param samples Character, specifying the name of samples to be compared.
#' @param main As in 'plot'.
#' @param ylab As in 'plot'.
#' @param filename Figure file name to create on disk. Default filename="NULL", which means
#' don't save the figure on disk.
#' @param width As in ggsave.
#' @param height As in ggsave.
#' @param ... Other available parameters in function 'ggsave'.
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
#' @seealso \code{\link{DensityView}}
#'
#' @examples
#' data(MLE_Data)
#' # Read beta score from gene summary table in MAGeCK MLE results
#' dd = ReadBeta(MLE_Data, organism="hsa")
#' ViolinView(dd, samples=c("D7_R1", "D7_R2", "PLX7_R1", "PLX7_R2"))
#' #or
#' ViolinView(dd[, 3:6])
#'
#'
#' @importFrom reshape melt
#' @importFrom ggsci scale_color_npg
#'
#' @export
#'

#===Distribution of beta scores======================================
ViolinView <- function(beta, samples=NULL, main=NULL,ylab="Beta Score",filename=NULL, width=5, height=4, ...){
  requireNamespace("reshape", quietly=TRUE) || stop("need reshape package")
  requireNamespace("ggsci", quietly=TRUE) || stop("need ggsci package")

  dd1=beta
  loginfo(paste("Violin plot for", main, ylab, "..."))
  if(!is.null(samples) && length(samples)>1){ dd1 = dd1[, samples]}

  dd1 = melt(dd1, id=NULL)
  #======
  p=ggplot(data=dd1,aes(x=variable,y=value,color=variable))
  p=p+geom_violin()+geom_boxplot(width=.1, outlier.colour=NA)
  #p=p+ylim(-1.5,1)
  p=p+scale_color_npg()
  p = p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank())
  p=p+theme(legend.position = "none")
  p=p+labs(x=NULL,y=ylab,title=main)

  if(!is.null(filename)){
    ggsave(plot=p,filename=filename, units = "in", dpi=600, width=width, height=height, ...)
  }
  return(p)
}
