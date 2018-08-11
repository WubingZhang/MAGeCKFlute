#' Density plot for gene beta scores in Control and Treatment
#'
#' Plot the density of gene beta scores in two samples.
#'
#' @docType methods
#' @name DensityView
#' @rdname DensityView
#'
#' @param beta Data frame, including \code{samples} as columns.
#' @param samples Character, specifying sample names in \code{beta}.
#' @param main As in 'plot'.
#' @param xlab As in 'plot'.
#' @param filename Figure file name to create on disk. Default filename="NULL", which means
#' don't save the figure on disk.
#' @param width As in ggsave.
#' @param height As in ggsave.
#' @param ... Other available parameters in ggsave.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#' @seealso \code{\link{ViolinView}}
#'
#' @examples
#' data(MLE_Data)
#' # Read beta score from gene summary table in MAGeCK MLE results
#' dd = ReadBeta(MLE_Data, organism="hsa")
#' DensityView(dd, samples=c("D7_R1", "D7_R2", "PLX7_R1", "PLX7_R2"))
#' #or
#' DensityView(dd[, c("D7_R1", "D7_R2", "PLX7_R1", "PLX7_R2")])
#'
#' @importFrom data.table melt
#' @importFrom ggsci scale_color_npg
#'
#' @export

#===Distribution of beta scores======================================
DensityView <- function(beta, samples=NULL, main=NULL,xlab="Beta Score",filename=NULL, width=5, height =4, ...){
  message(Sys.time(), " # Density plot for ", main, " ", xlab, " ...")
  if(!is.null(samples) && length(samples)>1){ beta = beta[, samples]}
  dd1 = data.table::melt(beta,id=NULL)
  if(!"variable" %in% colnames(dd1)){
    dd1$variable = colnames(beta)
  }
  #==========
  p=ggplot(data=dd1,aes(x=value,color=variable,group=variable))
  p=p+geom_density()
  # p=p+facet_wrap(~variable,nrow=1)
  p=p+scale_color_npg()
  p=p+labs(color=NULL)
  p=p+theme(legend.justification = c(1, 1), legend.position = c(0.99, 0.99))
  # p=p+theme(legend.text = element_text(size=8))
  p=p+labs(x=xlab, y="Density", title=main)
  p = p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(),
                legend.key = element_rect(fill = "transparent"))

  if(!is.null(filename)){
    ggsave(plot=p, filename=filename, units = "in", dpi=600, width=width, height=height, ...)
  }
  return(p)
}
