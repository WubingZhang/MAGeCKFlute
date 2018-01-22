#' Cluster and view cluster tree
#'
#' Cluster and view cluster tree
#'
#' @docType methods
#' @name hclustView
#' @rdname hclustView
#'
#' @param d A dissimilarity structure as produced by dist.
#' @param method The agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param label_cols A vector to be used as label's colors for the dendrogram.
#' @param bar_cols Either a vector or a matrix, which will be plotted as a colored bar.
#' @param main As in 'plot'.
#' @param xlab As in 'plot'.
#' @param horiz Logical indicating if the dendrogram should be drawn horizontally or not.
#' @param ... Arguments to be passed to methods, such as graphical parameters (see par).
#'
#' @return Just plot figure in device.
#'
#' @author Wubing Zhang
#'
#' @note See the vignette for an example of hclustView
#' Note that the source code of \code{hclustView} is very simple.
#' The source can be found by typing \code{MAGeCKFlute:::hclustView}
#' or \code{getMethod("hclustView")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/hclustView.R}
#' Users should find it easy to customize this function.
#'
#' @examples
#' label_cols = rownames(USArrests)
#' hclustView(dist(USArrests), label_cols=label_cols, bar_cols=label_cols)
#'
#' @importFrom dendextend labels_colors
#' @importFrom dendextend colored_bars
#'
#' @export
#'

hclustView <- function(d, method="average", label_cols=NULL, bar_cols=NULL, main=NA, xlab=NA, horiz = TRUE, ...){
  requireNamespace("dendextend")
  hc = hclust(dist(d), method = method)
  if(!is.null(label_cols)) label_cols = getCols(label_cols)
  if(!is.null(bar_cols)){
    if(is.null(dim(bar_cols))){
      bar_cols = getCols(bar_cols)
      mar4=1+max(nchar(hc$labels))/2
    }else{
      bar_cols = apply(bar_cols, 2, getCols)
      mar4=ncol(bar_cols)+max(nchar(hc$labels))/2
    }
  }else mar4 = max(nchar(hc$labels))/2

  par(mar=c(4, 2, 2, mar4), ...)  #
  dend=as.dendrogram(hc)
  if(!is.null(label_cols)) dendextend::labels_colors(dend) <- label_cols[order.dendrogram(dend)]
  plot(dend, main=main, xlab=xlab, horiz = horiz)
  if(!is.null(bar_cols)) dendextend::colored_bars(bar_cols, dend, horiz = horiz)
}
