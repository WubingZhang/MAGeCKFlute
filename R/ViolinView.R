#' Violin plot
#'
#' Plots the violin of beta scores in Control and Treatment samples.
#'
#' @docType methods
#' @name ViolinView
#' @rdname ViolinView
#' @aliases violinview
#'
#' @param beta Data frame, , including \code{samples} as columns.
#' @param samples Character, specifying the name of samples to be compared.
#' @param main As in 'plot'.
#' @param ylab As in 'plot'.
#' @param filename Figure file name to create on disk. Default filename="NULL", which means
#' don't save the figure on disk.
#' @param width As in ggsave.
#' @param height As in ggsave.
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
#' ViolinView(dd, samples=c("dmso", "plx"))
#' #or
#' ViolinView(dd[, c("dmso", "plx")])
#'
#'
#' @importFrom reshape2 melt
#'
#' @export
#'
ViolinView <- function(beta, samples=NULL, main=NULL, ylab="Beta Score",
                       filename=NULL, width=5, height=4, ...){

  if(!is.null(samples) && length(samples)>1){ beta = beta[, samples]}

  dd1 = reshape2::melt(beta, id.vars=NULL)
  if(!"variable" %in% colnames(dd1)){
    dd1$variable = colnames(beta)
  }
  #======
  p=ggplot(data=dd1,aes_string(x="variable",y="value",color="variable"))
  p=p+geom_violin()+geom_boxplot(width=.1, outlier.colour=NA)
  #p=p+ylim(-1.5,1)
  # p=p+scale_color_npg()
  p = p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(),
                legend.key = element_rect(fill = "transparent"))
  p=p+theme(legend.position = "none")
  p=p+labs(x=NULL,y=ylab,title=main)

  if(!is.null(filename)){
    ggsave(plot=p,filename=filename, units = "in", width=width, height=height, ...)
  }
  return(p)
}
