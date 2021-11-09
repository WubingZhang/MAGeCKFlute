#' Violin plot
#'
#' Violin plot showing the distribution of numeric vectors with the same length.
#'
#' @docType methods
#' @name ViolinView
#' @rdname ViolinView
#' @aliases violinview
#'
#' @param dat A data frame.
#' @param samples A character vector, specifying the columns in the \code{dat} for plotting.
#' @param main A character, specifying title.
#' @param ylab A character, specifying title of y-axis.
#' @param filename A character, specifying a file name to create on disk.
#' Set filename to be "NULL", if don't want to save the figure.
#' @param width Numeric, specifying width of figure.
#' @param height Numeric, specifying height of figure.
#' @param ... Other available parameters in function 'ggsave'.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#' @seealso \code{\link{DensityView}}
#'
#' @examples
#' file3 = file.path(system.file("extdata", package = "MAGeCKFlute"),
#' "testdata/mle.gene_summary.txt")
#' dd = ReadBeta(file3)
#' ViolinView(dd[, -1])
#'
#'
#' @importFrom reshape2 melt
#' @import ggplot2
#'
#' @export
#'
ViolinView <- function(dat, samples = NULL, main = NULL, ylab = "Score",
                       filename=NULL, width=5, height=4, ...){
  requireNamespace("ggplot2", quietly=TRUE) || stop("need ggplot2 package")
  requireNamespace("reshape2", quietly=TRUE) || stop("need reshape2 package")
  if(!is.null(samples) && length(samples)>1){
    dat = dat[, samples]
  }

  dd1 = reshape2::melt(dat, id.vars=NULL)
  if(!"variable" %in% colnames(dd1)){
    dd1$variable = colnames(dat)
  }
  ## Plotting
  p = ggplot(data=dd1, aes_string(x="variable",y="value",color="variable"))
  p = p + geom_violin()
  p = p + geom_boxplot(width=.1, outlier.colour=NA)
  p = p + theme_bw(base_size = 14)
  p = p + theme(plot.title = element_text(hjust = 0.5))
  p = p + theme(legend.position = "none")
  p = p + labs(x=NULL, y=ylab, title=main)
  if(!is.null(filename)){
    ggsave(plot=p,filename=filename, units = "in", width=width, height=height, ...)
  }
  return(p)
}
