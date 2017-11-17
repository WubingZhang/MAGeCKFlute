#' Density plot for gene beta scores in Control and Treatment
#'
#' Plot the density of gene beta scores in two samples.
#'
#' @docType methods
#' @name DensityView
#' @rdname DensityView
#'
#' @param beta data frame, which has columns of 'Gene', \code{ctrlname} and \code{treatname}.
#' @param ctrlname a character, specifying the name of control sample.
#' @param treatname a character, specifying the name of treatment sample.
#' @param main as in 'plot'.
#' @param xlab as in 'plot'.
#' @param filename figure file name to create on disk. Default filename="NULL", which means
#' don't save the figure on disk.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#' @note See the vignette for an example of density plot for beta score deviation.
#' Note that the source code of \code{DensityView} is very simple.
#' The source can be found by typing \code{MAGeCKFlute:::DensityView}
#' or \code{getMethod("DensityView")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/DensityView.R}
#' Users should find it easy to customize this function.
#'
#' @seealso \code{\link{DensityDiffView}}   \code{\link{KeggPathwayView}}
#' @seealso \code{\link{ViolinView}}   \code{\link{SquareView}}
#' @seealso \code{\link{CellCycleView}}  \code{\link{EnrichedView}}
#' @seealso \code{\link{EnrichedGSEView}}  \code{\link{MAView}}
#' @seealso \code{\link{RankView}}    \code{\link{ScatterView}}
#'
#' @examples
#' data(MLE_Data)
#' # Read beta score from gene summary table in MAGeCK MLE results
#' dd = ReadBeta(MLE_Data, ctrlName = "D7_R1", treatName = "PLX7_R1", organism="hsa")
#' #Density plot of all gene beta scores in control and treatment
#' DensityView(dd)
#'
#'
#' @importFrom reshape melt
#' @importFrom ggsci scale_color_npg
#'
#' @export

#===Distribution of beta scores======================================
DensityView <- function(beta, ctrlname="Control",treatname="Treatment",
                         main=NULL,xlab="Beta Score",filename=NULL){
  dd1 = beta
  loginfo(paste("Density plot for", main, xlab, "..."))
  dd1 = dd1[,c("Gene", ctrlname, treatname)]
  dd1 = melt(dd1,id="Gene")
  #==========
  p=ggplot(data=dd1,aes(x=value,color=variable,group=variable))
  p=p+geom_density()
  # p=p+facet_wrap(~variable,nrow=1)
  p=p+scale_color_npg()+theme_bw(14)+theme(plot.title = element_text(hjust = 0.5,size=12))
  p=p+labs(color=NULL)
  p=p+theme(legend.justification = c(0, 1), legend.position = c(0.01, 0.99))
  # p=p+theme(legend.text = element_text(size=8))
  p=p+labs(x=xlab, y="Density", title=main)

  if(!is.null(filename)){
    ggsave(plot=p,filename=filename,units = "in",width=300/100,height =270/100 )
  }
  return(p)
}
