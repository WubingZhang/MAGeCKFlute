#' MAplot of gene beta scores
#'
#' MAplot of gene beta scores in Control vs Treatment
#'
#' @docType methods
#' @name MAView
#' @rdname MAView
#'
#' @param beta Data frame, including \code{ctrlname} and \code{treatname} as columns.
#' @param ctrlname Character vector, specifying the name of control sample.
#' @param treatname Character vector, specifying the name of treatment sample.
#' @param main As in plot.
#' @param show.statistics Show statistics .
#' @param add.smooth Whether add a smooth line to the plot.
#' @param lty Line type for smooth line.
#' @param smooth.col Color of smooth line.
#' @param plot.method A string specifying the method to fit smooth line, which should be one of "loess" (default), "lm", "glm" and "gam".
#' @param filename Figure file name to create on disk. Default filename="NULL", which means
#' don't save the figure on disk.
#' @param width As in ggsave.
#' @param height As in ggsave.
#' @param ... Other available parameters in function 'ggsave'.
#'
#' @author Wubing Zhang
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#'
#' @examples
#' data(mle.gene_summary)
#' # Read beta score from gene summary table in MAGeCK MLE results
#' dd = ReadBeta(mle.gene_summary, organism="hsa")
#' MAView(dd, ctrlname = "dmso", treatname = "plx")
#'
#' @export

MAView <- function(beta, ctrlname="Control",treatname="Treatment", main=NULL,
                    show.statistics = TRUE, add.smooth = TRUE, lty = 1, smooth.col = "red",
                    plot.method = c("loess", "lm", "glm", "gam"),
                    filename=NULL, width=5, height=4, ...){
  dd = beta
  dd[is.na(dd)] = 0
  message(Sys.time(), " # MAplot for ", main, " beta scores ...")
  A = rowMeans(dd[,c(ctrlname, treatname)])
  M = rowMeans(dd[,treatname,drop= FALSE])-rowMeans(dd[,ctrlname,drop= FALSE])
  subset = sample(1:length(M), min(c(10000, length(M))))
  A = A[subset]
  M = M[subset]
  Mid = paste("Median: ", round(median(M),4), sep="")
  IQR = round(quantile(M, 0.75) - quantile(M, 0.25), 4)
  IQR = paste("IQR: ", IQR, sep="")
  gg = data.frame(M=M, A=A)
  p = ggplot(gg, aes(x=A, y=M))
  p = p + geom_point(shape=1, size=0.5, alpha = 0.6)
  p = p + geom_hline(yintercept = 0, color="blue")
  if(add.smooth)
    p = p + geom_smooth(method = plot.method[1], color=smooth.col, linetype=lty)
  p = p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(),
                legend.key = element_rect(fill = "transparent"))
  p = p + labs(title=main)
  if(show.statistics){
    xmax = max(gg$A)
    ymax = max(gg$M)
    p = p + annotate("text", color="black", x=xmax, y=ymax, hjust = 1, vjust = 1,
                     label=paste(Mid, IQR, sep="\n"))
  }
  if(!is.null(filename)){
    ggsave(plot=p, filename=filename, units = "in", width=width, height =height, ...)
  }
  return(p)
}
