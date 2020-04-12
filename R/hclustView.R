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
#' @return Plot figure on open device.
#'
#' @author Wubing Zhang
#'
#' @examples
#' label_cols = rownames(USArrests)
#' hclustView(dist(USArrests), label_cols=label_cols, bar_cols=label_cols)
#'
#' @importFrom dendextend labels_colors
#' @importFrom dendextend colored_bars
#' @importFrom graphics par plot
#' @export
#'

hclustView <- function(d, method="average", label_cols=NULL, bar_cols=NULL, main=NA, xlab=NA, horiz = TRUE, ...){
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

  graphics::par(mar=c(4, 2, 2, mar4), ...)  #
  dend=as.dendrogram(hc)
  if(!is.null(label_cols)) dendextend::labels_colors(dend) <- label_cols[order.dendrogram(dend)]
  graphics::plot(dend, main=main, xlab=xlab, horiz = horiz)
  if(!is.null(bar_cols)) dendextend::colored_bars(bar_cols, dend, horiz = horiz)
}

#' Map values to colors
#'
#' @docType methods
#' @name getCols
#' @rdname getCols
#'
#' @param x A numeric vector.
#' @param palette diverge, rainbow, sequential
#'
#' @return A vector of colors corresponding to input vector.
#'
#' @author Wubing Zhang
#'
#'
#' @importFrom scales gradient_n_pal
#' @importFrom scales brewer_pal
#'
#' @examples
#' getCols(1:4)
#' @export
#'
getCols <- function(x, palette=1){
  cols = scales::gradient_n_pal(scales::brewer_pal("qual", palette)(8))(seq(0, 1, length.out = length(unique(x))))
  names(cols) = unique(x)
  return(cols[x])
}
