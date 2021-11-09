#' Density plot
#'
#' Plot the distribution of numeric vectors with the same length.
#'
#' @docType methods
#' @name DensityView
#' @rdname DensityView
#'
#' @param dat A data frame.
#' @param samples A character vector, specifying columns in \code{dat} for plotting.
#' @param main A character, specifying title.
#' @param xlab A character, specifying title of x-axis.
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
#' @seealso \code{\link{ViolinView}}
#'
#' @examples
#' file3 = file.path(system.file("extdata", package = "MAGeCKFlute"),
#' "testdata/mle.gene_summary.txt")
#' dd = ReadBeta(file3)
#' DensityView(dd, samples=c("Pmel1_Ctrl", "Pmel1"))
#' #or
#' DensityView(dd[,-1])
#'
#' @importFrom reshape2 melt
#'
#' @export
DensityView <- function(dat, samples = NULL,
                        main = NULL, xlab = "Score",
                        filename = NULL, width = 5, height = 4, ...){
  if(!is.null(samples) && length(samples)>0){
    dat = dat[, samples, drop = FALSE]
  }
  dd1 = reshape2::melt(dat,id=NULL)
  if(!"variable" %in% colnames(dd1)){
    dd1$variable = colnames(dat)
  }
  #==========
  p = ggplot(data=dd1, aes_string(x="value", color="variable", group="variable"))
  p = p + geom_density()
  p = p + labs(color=NULL)
  p = p + theme(legend.justification = c(1, 1), legend.position = c(0.99, 0.99))
  # p=p+theme(legend.text = element_text(size=8))
  p = p + labs(x=xlab, y="Density", title=main)
  p = p + theme_bw(base_size = 14)
  p = p + theme(plot.title = element_text(hjust = 0.5))

  if(!is.null(filename)){
    ggsave(plot=p, filename=filename, units = "in",
           width=width, height=height, ...)
  }
  return(p)
}



#' Density plot
#'
#' Plot the distribution of score differences between treatment and control.
#'
#' @docType methods
#' @name DensityDiffView
#'
#' @param dat A data frame.
#' @param ctrlname A character, specifying the control samples.
#' @param treatname A character, specifying the treatment samples.
#' @param main A character, specifying title.
#' @param filename A character, specifying a file name to create on disk.
#' Set filename to be "NULL", if don't want to save the figure.
#' @param width Numeric, specifying width of figure.
#' @param height Numeric, specifying height of figure.
#' @param ... Other parameters in ggsave.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#'
#' @examples
#' file3 = file.path(system.file("extdata", package = "MAGeCKFlute"),
#' "testdata/mle.gene_summary.txt")
#' dd = ReadBeta(file3)
#' # Density plot of beta score deviation between control and treatment
#' DensityDiffView(dd, ctrlname = "Pmel1_Ctrl", treatname = "Pmel1")
#'
#'
#' @export

DensityDiffView <- function(dat, ctrlname="Control", treatname="Treatment", main=NULL,
                            filename=NULL, width = 5, height = 4, ...){
  d = dat
  d$Diff = rowMeans(d[,treatname,drop=FALSE])-rowMeans(d[,ctrlname,drop=FALSE])
  d$r = rnorm(length(d$Diff), mean=0, sd=sd(d$Diff)-0.01)
  p = ggplot(d, aes_string("Diff"))
  p = p + geom_density(colour="black")
  # p = p + geom_histogram(aes(y = ..density..), fill="gray90", binwidth=0.02)
  p = p + geom_density(aes_string("r"), linetype="dashed", colour="red")
  p = p + geom_vline(xintercept = 0,linetype="dashed")
  p = p + theme_bw(base_size = 14)
  p = p + theme(plot.title = element_text(hjust = 0.5))
  p = p + labs(x="Treatment vs control", y="Density", title=main)

  if(!is.null(filename)){
    ggsave(plot=p, filename=filename, units = "in", width=width, height=height, ...)
  }
  return(p)
}
