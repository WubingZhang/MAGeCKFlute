#' Calculate the similarity between samples and plot heatmap
#'
#' Calculate the similarity between samples and plot heatmap
#'
#' @docType methods
#' @name HeatmapView
#' @rdname HeatmapView
#' @aliases heatmapview
#'
#' @param beta Data frame or matrix, in which each column represents one sample.
#' @param method Character, One of "pearson", "kendall", "spearman", "euclidean", "maximum",
#' "manhattan", "canberra", "binary", or "minkowski".
#' @param breaks The same as that in pheatmap
#' @param cluster_rows The same as that in pheatmap
#' @param cluster_cols The same as that in pheatmap
#' @param legend The same as that in pheatmap
#' @param main The same as that in pheatmap
#' @param fontsize The same as that in pheatmap
#' @param display_numbers The same as that in pheatmap
#' @param filename The same as that in pheatmap
#' @param width The same as that in pheatmap
#' @param height The same as that in pheatmap
#' @param ... Other parameters in pheatmap
#'
#' @return The same as pheatmap
#'
#' @author Wubing Zhang
#'
#' @seealso \code{\link[pheatmap]{pheatmap}}
#'
#' @examples
#' data(MLE_Data)
#' dd = ReadBeta(MLE_Data, organism="hsa")
#' dd = dd[,3:ncol(dd)]
#' HeatmapView(dd, method = "pearson")
#'
#' @import pheatmap
#' @importFrom grDevices colorRampPalette
#'
#' @export
#'

HeatmapView <- function(beta, method = "pearson", breaks = NA, cluster_rows = TRUE, cluster_cols = TRUE,
                        legend = TRUE, main = NA, fontsize = 10, display_numbers = TRUE,
                        filename = NA, width = NA, height = NA, ...){
  if(method%in%c("pearson", "kendall", "spearman")){
    mat = cor(beta, method=method)
    if(is.na(breaks[1])) breaks = seq(-1, 1, length.out = 200)
    color = rev(colorRampPalette(RColorBrewer::brewer.pal(7, "RdBu"))(length(breaks)-1))
  }else if(method %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")){
    mat = as.matrix(dist(t(beta), method=method))
    if(is.na(breaks[1])) breaks = seq(min(mat),max(mat), length.out = 200)
    color = rev(colorRampPalette(RColorBrewer::brewer.pal(7, "Reds"))(length(breaks)-1))
  }
  pheatmap::pheatmap(mat, color=color, breaks=breaks, cluster_rows=cluster_rows,
           cluster_cols=cluster_cols, legend=legend, main = main, fontsize=fontsize, fontfamily = "Helvetica",
           display_numbers=display_numbers, filename=filename, width=width, height = height,...)
}
