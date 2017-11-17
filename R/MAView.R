#' MAplot of gene beta scores
#'
#' MAplot of gene beta scores in Control vs Treatment
#'
#' @docType methods
#' @name MAView
#' @rdname MAView
#'
#' @param beta data frame, which has columns of 'Gene', \code{ctrlname} and \code{treatname}.
#' @param ctrlname character vector, specifying the name of control sample.
#' @param treatname character vector, specifying the name of treatment sample.
#' @param main as in "ma.plot".
#' @param subset as in "ma.plot".
#' @param show.statistics as in "ma.plot".
#' @param span as in "ma.plot".
#' @param family.loess as in "ma.plot".
#' @param cex as in "ma.plot".
#' @param cex.lab as in "ma.plot".
#' @param cex.axis as in "ma.plot".
#' @param cex.main as in "ma.plot".
#' @param plot.method as in "ma.plot".
#' @param add.loess as in "ma.plot".
#' @param lwd as in "ma.plot".
#' @param lty as in "ma.plot".
#' @param loess.col as in "ma.plot".
#' @param filename figure file name to create on disk. Default filename="NULL", which means
#' no output.
#'
#' @author Wubing Zhang
#'
#' @note See the vignette for an example of MAView.
#' Note that the source code of \code{MAView} is very simple.
#' The source can be found by typing \code{MAGeCKFlute:::MAView}
#' or \code{getMethod("MAView")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/MAView.R}
#' Users should find it easy to customize this function.
#'
#' @seealso \code{\link{DensityDiffView}}   \code{\link{DensityView}}
#' @seealso \code{\link{ViolinView}}   \code{\link{SquareView}}
#' @seealso \code{\link{CellCycleView}}  \code{\link{EnrichedView}}
#' @seealso \code{\link{EnrichedGSEView}}  \code{\link{KeggPathwayView}}
#' @seealso \code{\link{RankView}}    \code{\link{ScatterView}}
#'
#' @examples
#' data(MLE_Data)
#' # Read beta score from gene summary table in MAGeCK MLE results
#' dd = ReadBeta(MLE_Data, ctrlName = "D7_R1", treatName = "PLX7_R1", organism="hsa")
#' MAView(dd, cex=1)
#'
#' @importFrom affy ma.plot
#' @importFrom graphics par plot.new
#' @importFrom grid viewport
#'
#' @export

MAView <- function(beta, ctrlname="Control",treatname="Treatment", main="Negative control normalized", subset = sample(1:length(M), min(c(10000, length(M)))),
                    show.statistics = TRUE, span = 2/3, family.loess = "gaussian",
                    cex = 1, cex.lab=1.2, cex.axis=1, cex.main=1.2, plot.method = c("normal","smoothScatter","add"),
                    add.loess = TRUE, lwd = 1, lty = 1, loess.col = "red",filename=NULL){
  requireNamespace("graphics", quietly=TRUE) || stop("need graphics package")
  requireNamespace("grid", quietly=TRUE) || stop("need grid package")
  requireNamespace("affy", quietly=TRUE) || stop("need affy package")

  dd=beta
  loginfo(paste("MAplot for", main, "beta scores ..."))
  A = rowMeans(dd[,c(ctrlname, treatname)])
  M = rowMeans(dd[,treatname,drop= FALSE])-rowMeans(dd[,ctrlname,drop= FALSE])
  par(cex.axis=cex.axis, cex.lab=cex.lab,cex.main=cex.main)
  ma.plot(A, M, main=main, subset=subset, show.statistics=show.statistics, span=span,
          family.loess=family.loess, cex=cex, plot.method=plot.method, add.loess=add.loess,
          lwd=lwd, lty=lty, loess.col=loess.col)
  if(!is.null(filename)){
    png(filename,units = "in",width=500/100,height =450/100, res=300)
    par(cex.axis=cex.axis, cex.lab=cex.lab,cex.main=cex.main)
    ma.plot(A, M, main=main, subset=subset, show.statistics=show.statistics, span=span,
            family.loess=family.loess, cex=cex, plot.method=plot.method, add.loess=add.loess,
            lwd=lwd, lty=lty, loess.col=loess.col)
    dev.off()
  }
}
