#' Estimate cell cycle time for all samples compared to control sample and view.
#'
#' Estimate cell cycle time in different samples by linear fitting of beta scores, and plot fitting lines,
#' in which x-axis is control beta score and y-axis is beta score of all samples.
#'
#' @docType methods
#' @name CellCycleView
#' @rdname CellCycleView
#' @aliases CellCycle,MAGeCKFlute-method
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
#' data(mle.gene_summary)
#' # Read beta score from gene summary table in MAGeCK MLE results
#' dd = ReadBeta(mle.gene_summary, organism="hsa")
#' CellCycleView(dd, ctrlname = "dmso", treatname = "plx")
#'
#' @importFrom data.table melt
#' @importFrom ggsci scale_color_npg
#' @export

#===Distribution of beta scores======================================
CellCycleView <- function(beta, ctrlname, treatname, main=NULL, filename=NULL, width=5, height = 4, ...){
  message(Sys.time(), " # Cell cycle fitting of treatments compaired to ",
                paste(ctrlname, collapse = "&"))

  requireNamespace("data.table", quietly=TRUE) || stop("need data.table package")
  dd2 = data.frame(x = rowMeans(beta[,ctrlname,drop=FALSE]), y = rowMeans(beta[,treatname,drop=FALSE]))
  p = ggplot(dd2, aes(x, y))
  p = p + geom_jitter(size = 0.1, alpha=0.8, color="#1f78b4")
  p = p + geom_abline(slope = 1, intercept = 0, color="#47484a", linetype=2, size=0.3)
  p = p + geom_smooth(method='lm', se=TRUE, size=0.5, color="#e41a1c")
  p = p + labs(x="Control", y="Treatment", title=main, color=NULL)
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


