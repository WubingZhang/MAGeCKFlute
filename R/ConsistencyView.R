#' Visualize the estimate cell cycle compared to control.
#'
#' Estimate cell cycle time in different samples by linear fitting of beta scores.
#'
#' @docType methods
#' @name ConsistencyView
#' @rdname ConsistencyView
#'
#' @param beta Data frame, which has columns of \code{ctrlname} and other samples.
#' @param ctrlname A character, specifying the names of control samples.
#' @param treatname A character, specifying the name of treatment samples.
#' @param main As in 'plot'.
#' @param filename Figure file name to create on disk. Default filename="NULL", which means
#' no output.
#' @param width As in ggsave.
#' @param height As in ggsave.
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
#' ConsistencyView(dd, ctrlname = "dmso", treatname = "plx")
#'
#' @export

ConsistencyView <- function(beta, ctrlname, treatname, main=NULL,
                          filename=NULL, width=5, height = 4, ...){
  dd2 = data.frame(x = rowMeans(beta[,ctrlname,drop=FALSE]),
                   y = rowMeans(beta[,treatname,drop=FALSE]))
  p = ScatterView(dd2, "x", "y", color="#1f78b4")
  p = p + geom_abline(slope = 1, intercept = 0, color="gray50", linetype=2, size=0.8)
  p = p + geom_smooth(method='lm', formula = y ~ x, se=TRUE, size=0.5, color="#e41a1c")
  p = p + labs(x = ctrlname[1], y = treatname[1], title=main, color=NULL)
  p = p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(),
                legend.key = element_rect(fill = "transparent"))
  p = p + theme(legend.position = "none")
  if(!is.null(filename)){
    ggsave(plot=p, filename=filename, units = "in", width=width, height = height, ...)
  }
  return(p)
}


