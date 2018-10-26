#' Visualize the correlation between two object
#'
#' @docType methods
#' @name CorrView
#' @rdname CorrView
#'
#' @param gg A data frame.
#' @param x A character, indicating column (in countSummary) of x-axis.
#' @param y A character, indicating column (in countSummary) of y-axis.
#' @param smoothMethod A character, indicating fill color of all bars.
#' @param main A charater, specifying the figure title.
#' @param xlab A character, specifying the title of x-axis.
#' @param ylab, A character, specifying the title of y-axis.
#' @param filename Figure file name to create on disk. Default filename="NULL", which means
#' don't save the figure on disk.
#' @param width As in ggsave.
#' @param height As in ggsave.
#' @param ... Other available parameters in ggsave.
#'
#' @return An object created by \code{ggplot}, which can be assigned
#' and further customized.
#'
#' @examples
#' gg = data.frame(x = rnorm(50), y = rnorm(50))
#' CorrView(gg, x="x", y="y")
#'

#' @import ggplot2
#' @export

CorrView <- function(gg, x, y, smoothMethod = "lm",
                     main = NULL, xlab = NULL, ylab = NULL,
                     filename = NULL, width = 5, height = 4, ...){
  gg$x = gg[, x]
  gg$y = gg[, y]

  gg = gg[is.finite(gg$x) & is.finite(gg$y), ]
  tmp = cor.test(gg$x, gg$y)
  corr = paste('Correlation = ', round(tmp$estimate,3),
               "\nP.value = ", signif(tmp$p.value,3))

  p = ggplot(gg, aes(x, y))
  p = p + geom_point(size=0.5, color="blue")
  p = p + geom_smooth(method=smoothMethod, se=FALSE, size=0.5, color = "#c12603")
  # p = p + geom_abline(slope = 1, intercept = 0, color="gray50", linetype=2)
  p = p + labs(x=xlab, y=ylab, title = main)
  p = p + theme(legend.key = element_rect(fill = "transparent"))
  p = p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank())
  p = p + theme(legend.justification = c(0, 1), legend.position = c(0.01, 0.99))
  if(tmp$estimate>0){
    p = p + annotate("text", x = max(gg$x), y = min(gg$y),
                     label = corr, hjust=1, vjust=0)

  }else{
    p = p + annotate("text", x = max(gg$x), y = max(gg$y),
                     label = corr, hjust=1, vjust=1)
  }
  if(!is.null(filename)){
    ggsave(plot=p, filename=filename, units = "in", width=width, height=height, ...)
  }
  return(p)
}
