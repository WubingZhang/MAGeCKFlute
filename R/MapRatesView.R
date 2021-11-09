#' View mapping ratio
#'
#' View mapping ratio of each sample
#'
#' @docType methods
#' @name MapRatesView
#' @rdname MapRatesView
#'
#' @param countSummary A data frame, which contains columns of
#' `Label`, `Reads`, and `Mapped`
#' @param Label A character, indicating column (in countSummary) of sample names.
#' @param Reads A character, indicating column (in countSummary) of total reads.
#' @param Mapped A character, indicating column (in countSummary) of mapped reads.
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
#' MapRatesView(countsummary)
#'
#' @import ggplot2
#' @export

MapRatesView <- function(countSummary,
                         Label = "Label",
                         Reads = "Reads",
                         Mapped = "Mapped",
                         filename = NULL,
                         width = 5, height = 4,
                         ...){
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package \"scales\" is required. Please install it.", call. = FALSE)
  }
  gg = data.frame(Label=rep(countSummary[, Label], 2),
                  read=rep(countSummary[, Reads],2),
                  count=c(countSummary[, Mapped],
                          countSummary[, Reads]-countSummary[, Mapped]),
                  category=factor(rep(c("mapped", "unmapped"),
                                      c(nrow(countSummary), nrow(countSummary))),
                                  levels = c("unmapped", "mapped")))
  gg$percent = paste0(round(gg$count*100/gg$read, 1), "%")
  gg$pos = ceiling(gg$count/2)
  gg$pos[gg$category=="unmapped"] = ceiling(gg$pos[gg$category=="unmapped"] +
                                              gg$pos[gg$category=="mapped"]*2)

  fill = c("#9BC7E9", "#1C6DAB")
  p <- ggplot(gg)
  p = p + geom_bar(aes_string(y = "count", x = "Label", fill = "category"),
                   stat="identity", width=0.8, alpha=0.9)
  p = p + geom_text(aes_string(x = "Label", y = "pos", label = "percent"), size=4)
  p = p + labs(x=NULL, y="Reads", title="Mapping ratio")
  p = p + scale_y_continuous(expand = c(0,0))
  p = p + scale_fill_manual(values=fill)
  p = p + theme(legend.title = element_blank())
  p = p + theme_bw(base_size = 14)
  p = p + theme(plot.title = element_text(hjust = 0.5))
  if(!is.null(filename)){
    ggsave(plot=p, filename=filename, units = "in", width=width, height=height, ...)
  }
  return(p)
}
