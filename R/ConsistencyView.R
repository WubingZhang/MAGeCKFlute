#' Visualize the estimate cell cycle compared to control.
#'
#' Estimate cell cycle time in different samples by linear fitting of beta scores.
#'
#' @docType methods
#' @name ConsistencyView
#' @rdname ConsistencyView
#'
#' @param dat A data frame.
#' @param ctrlname A character, specifying the names of control samples.
#' @param treatname A character, specifying the names of treatment samples.
#' @param main A character, specifying title.
#' @param filename A character, specifying a file name to create on disk.
#' Set filename to be "NULL", if don't want to save the figure.
#' @param width Numeric, specifying width of figure.
#' @param height Numeric, specifying height of figure.
#' @param ... Other available parameters in ggsave.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#' @examples
#' file3 = file.path(system.file("extdata", package = "MAGeCKFlute"),
#' "testdata/mle.gene_summary.txt")
#' dd = ReadBeta(file3)
#' ConsistencyView(dd, ctrlname = "Pmel1_Ctrl", treatname = "Pmel1")
#'
#' @export

ConsistencyView <- function(dat, ctrlname, treatname, main=NULL,
                          filename=NULL, width=5, height = 4, ...){
  dd2 = data.frame(x = rowMeans(dat[,ctrlname,drop=FALSE]),
                   y = rowMeans(dat[,treatname,drop=FALSE]))
  p = ScatterView(dd2, "x", "y", color="#1f78b4")
  p = p + geom_abline(slope = 1, intercept = 0, color="gray50", linetype=2, size=0.8)
  p = p + geom_smooth(method='lm', formula = y ~ x, se=TRUE, size=0.5, color="#e41a1c")
  p = p + labs(x = ctrlname[1], y = treatname[1], title=main, color=NULL)
  p = p + theme_bw(base_size = 14)
  p = p + theme(plot.title = element_text(hjust = 0.5))
  p = p + theme(legend.position = "none")
  if(!is.null(filename)){
    ggsave(plot=p, filename=filename, units = "in", width=width, height = height, ...)
  }
  return(p)
}


