#' Density plot for beta score deviation between Control and Treatment
#'
#' Plot the density of beta score deviation between two samples.
#'
#' @docType methods
#' @name DensityDiffView
#'
#' @param beta Data frame, including \code{ctrlname} and \code{treatname} as columns.
#' @param ctrlname A character, specifying the name of control sample.
#' @param treatname A character, specifying the name of treatment sample.
#' @param main As in 'plot'.
#' @param filename Figure file name to create on disk. Default filename="NULL", which means no output.
#' @param width As in ggsave.
#' @param height As in ggsave.
#' @param ... Other parameters in ggsave.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#'
#' @examples
#' data(mle.gene_summary)
#' # Read beta score from gene summary table in MAGeCK MLE results
#' dd = ReadBeta(mle.gene_summary, organism="hsa")
#' # Density plot of beta score deviation between control and treatment
#' DensityDiffView(dd, ctrlname = "dmso", treatname = "plx")
#'
#'
#' @export

#===Distribution of beta scores======================================

DensityDiffView <- function(beta, ctrlname="Control", treatname="Treatment", main=NULL,
                            filename=NULL, width = 5, height = 4, ...){

  message(Sys.time(), " # Density plot for ", main, " treat-control beta scores...")
  d=beta
  d$Diff=rowMeans(d[,treatname,drop=FALSE])-rowMeans(d[,ctrlname,drop=FALSE])
  d$r <- rnorm(length(d$Diff), mean=0, sd=sd(d$Diff)-0.01)
  p=ggplot(d,aes(x=Diff))
  p=p+geom_histogram(aes(y = ..density..),fill="gray90",binwidth=0.02)
  p=p+geom_density(colour="black")
  p=p+geom_density(aes(x=r,y=..density..),linetype="dashed",colour="red")
  p=p+geom_vline(xintercept = 0,linetype="dashed")
  p = p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(),
                legend.key = element_rect(fill = "transparent"))
  p=p+labs(x="Treat-Control Beta Score",y="Density",title=main)
  #+ggtitle("Normalization with")

  if(!is.null(filename)){
    ggsave(plot=p, filename=filename, units = "in", width=width, height=height, ...)
  }
  return(p)
}
