#' Identical bar plot
#'
#' Identical bar plot
#'
#' @docType methods
#' @name IdentBarView
#' @rdname IdentBarView
#'
#' @param gg A data frame.
#' @param x A character, indicating column (in countSummary) of x-axis.
#' @param y A character, indicating column (in countSummary) of y-axis.
#' @param fill A character, indicating fill color of all bars.
#' @param main A charater, specifying the figure title.
#' @param xlab A character, specifying the title of x-axis.
#' @param ylab, A character, specifying the title of y-axis.
#' @param filename Figure file name to create on disk. Default filename="NULL", which means
#' don't save the figure on disk.
#' @param width As in ggsave.
#' @param height As in ggsave.
#' @param ... Other available parameters in ggsave.
#'
#' @author Wubing Zhang
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#'
#' @examples
#' file4 = file.path(system.file("extdata", package = "MAGeCKFlute"),
#' "testdata/countsummary.txt")
#' countsummary = read.delim(file4, check.names = FALSE)
#' IdentBarView(countsummary, x="Label", y="Reads")
#'
#' @import ggplot2
#' @export

IdentBarView <- function(gg, x = "x", y = "y", fill = c("#CF3C2B", "#394E80"),
                         main = NULL, xlab = NULL, ylab = NULL,
                         filename = NULL, width = 5, height = 4, ...){
  gg$x = gg[, x]
  gg$y = gg[, y]
  p <- ggplot(gg)
  p = p + geom_bar(aes(x, y), stat="identity", width=0.6, fill = fill[1], alpha=0.8)
  p = p + labs(x=xlab, y=ylab, title=main)
  p = p + scale_y_continuous(expand = c(0,0))
  p = p + theme_bw(base_size = 14)
  p = p + theme(plot.title = element_text(hjust = 0.5, size=18),
                axis.text.x=element_text(angle = 45, hjust=1, vjust = 1))

  if(!is.null(filename)){
    ggsave(plot=p, filename=filename, units = "in", width=width, height=height, ...)
  }
  return(p)
}
