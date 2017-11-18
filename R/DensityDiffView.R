#' Density plot for beta score deviation between Control and Treatment
#'
#' Plot the density of beta score deviation between two samples.
#'
#' @docType methods
#' @name DensityDiffView
#'
#' @param beta data frame, which has columns of 'Gene', \code{ctrlname} and \code{treatname}.
#' @param ctrlname a character, specifying the name of control sample.
#' @param treatname a character, specifying the name of treatment sample.
#' @param main as in 'plot'.
#' @param filename figure file name to create on disk. Default filename="NULL", which means
#' no output.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#' @note See the vignette for an example of density plot for beta score deviation.
#' Note that the source code of \code{DensityDiffView} is very simple.
#' The source can be found by typing \code{MAGeCKFlute:::DensityDiffView}
#' or \code{getMethod("DensityDiffView")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/DensityDiffView.R}
#' Users should find it easy to customize this function.
#'
#' @seealso \code{\link{KeggPathwayView}}   \code{\link{DensityView}}
#' @seealso \code{\link{ViolinView}}   \code{\link{SquareView}}
#' @seealso \code{\link{CellCycleView}}  \code{\link{EnrichedView}}
#' @seealso \code{\link{EnrichedGSEView}}  \code{\link{MAView}}
#' @seealso \code{\link{RankView}}    \code{\link{ScatterView}}
#'
#' @examples
#' data(MLE_Data)
#' # Read beta score from gene summary table in MAGeCK MLE results
#' dd = ReadBeta(MLE_Data, ctrlName = "D7_R1", treatName = "PLX7_R1", organism="hsa")
#' # Density plot of beta score deviation between control and treatment
#' DensityDiffView(dd)
#'
#'
#' @export

#===Distribution of beta scores======================================

DensityDiffView <- function(beta, ctrlname="Control",
                         treatname="Treatment",main=NULL,filename=NULL){

  loginfo(paste("Density plot for", main, "treat-control beta scores..."))
  d=beta
  d$Diff=rowMeans(d[,treatname,drop=FALSE])-rowMeans(d[,ctrlname,drop=FALSE])
  d$r <- rnorm(length(d$Diff), mean=0, sd=sd(d$Diff)-0.01)
  p=ggplot(d,aes(x=Diff))
  p=p+geom_histogram(aes(y = ..density..),fill="gray90",binwidth=0.02)
  p=p+geom_density(colour="black")
  p=p+geom_density(aes(x=r,y=..density..),linetype="dashed",colour="red")
  p=p+geom_vline(xintercept = 0,linetype="dashed")
  p=p+theme_bw(14)+theme(plot.title = element_text(hjust = 0.5,size=12))
  p=p+labs(x="Treat-Control Beta Score",y="Density",title=main)
  #+ggtitle("Normalization with")

  if(!is.null(filename)){
    ggsave(plot=p,filename=filename,units = "in",width=300/100,
           height =270/100 )
  }
  return(p)
}
