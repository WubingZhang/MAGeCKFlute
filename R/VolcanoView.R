#' Volcano View
#'
#' Volcano plot for differential analysis.
#'
#' @docType methods
#' @name VolcanoView
#' @rdname VolcanoView
#'
#' @param df A data frame.
#' @param x A character, specifying the x-axis in Volcanno figure, 'logFC' (default).
#' @param y A character, specifying the y-axis in Volcanno figure, 'adj.P.Val' (default).
#' log10 transformation will be done automatically.
#' @param Label A character, specifying dots to be labeled on the figure.
#' @param top An integer, specifying the number of top significant genes to be labeled.
#' @param topnames A character vector, indicating positive/negative controls to be labeled.
#' @param x_cutoff Numeric, specifying cutoff of the x-axis.
#' @param y_cutoff Numeric, specifying cutoff of the y-axis.
#' @param mycolour A color vector, specifying colors of non-significant, significantly up and down-regulated genes.
#' @param alpha Numeric, parameter in ggplot.
#' @param force Numeric, Parameter for geom_text_repel. Force of repulsion between overlapping text labels.
#' @param main A character, specifying title.
#' @param xlab A character, specifying title of x-axis.
#' @param ylab A character, specifying title of y-axis.
#' @param filename A character, specifying a file name to create on disk. Set filename to be "NULL",
#' if don't want to save the figure.
#' @param width Numeric, specifying width of figure.
#' @param height Numeric, specifying height of figure.
#' @param ... Other available parameters in ggsave.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#' @examples
#' file1 = file.path(system.file("extdata", package = "MAGeCKFlute"),
#' "testdata/rra.gene_summary.txt")
#' gdata = ReadRRA(file1)
#' VolcanoView(gdata, x = "Score", y = "FDR", Label = "id")
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @export

VolcanoView <- function(df, x = "logFC", y = "adj.P.Val",
                        Label = NA, top = 5, topnames = NULL,
                        x_cutoff = log2(1.5), y_cutoff = 0.05,
                        mycolour = c("gray80", "#e41a1c", "#377eb8"),
                        alpha = 0.6, force = 0.1, main = NULL,
                        xlab = "log2FC", ylab = "-log10(FDR)",
                        filename = NULL, width = 4, height = 2.5,
                        ...){
  requireNamespace("ggplot2", quietly=TRUE) || stop("need ggplot2 package")
  requireNamespace("ggrepel", quietly=TRUE) || stop("need ggrepel package")
    gg = df[, c(x, y)]
    gg$group="no"
    gg$group[gg[,x]>x_cutoff & gg[,y]<y_cutoff] = "up"
    gg$group[gg[,x]< -x_cutoff & gg[,y]<y_cutoff] = "down"

    gg[, y] = -log10(gg[, y])
    if(!(top==0 & is.null(topnames))){
      gg$Label = rownames(gg)
      if(!is.na(Label)) gg$Label = df[, Label]
      gg = gg[order(gg[,y], abs(gg[,x]), decreasing = TRUE), ]
      idx1 = idx2 = c()
      if(top>0){
        idx1 = which(gg$group=="up")[1:min(top, sum(gg$group=="up"))]
        idx2 = which(gg$group=="down")[1:min(top, sum(gg$group=="down"))]
      }
      idx = unique(c(idx1, idx2, which(gg$Label %in% topnames)))
      gg$Label = as.character(gg$Label)
      gg$Label[setdiff(1:nrow(gg), idx)] = ""
      # gg$Label = factor(gg$Label, levels = setdiff(unique(gg$Label), ""))
    }
    gg$color = gg$group
    gg$color[gg$Label!=""] = "black"
    mycolour = c(mycolour, "black")
    names(mycolour) = c("no", "up", "down", "black")
    #=========
    p = ggplot(gg, aes_string(x=x, y=y, label="Label"))
    p = p + geom_point(aes_string(fill="group"), shape = 21, alpha=alpha)
    p = p + geom_point(aes_string(colour="color"), shape = 21, alpha=alpha)
    p = p + scale_color_manual(values=mycolour)
    p = p + scale_fill_manual(values=mycolour)
    p = p + geom_hline(yintercept = -log10(y_cutoff), linetype = "dotted")
    p = p + geom_vline(xintercept = c(-x_cutoff, x_cutoff), linetype = "dotted")
    p = p + xlim(min(gg[,x])-0.001, max(gg[,x])+0.001)
    p = p + labs(x=xlab, y=ylab, title=main)
    p = p + theme_bw(base_size = 14)
    p = p + theme(plot.title = element_text(hjust = 0.5))
    p = p + theme(legend.position = "none")
    if(!(top==0 & is.null(topnames))){
      p = p + ggrepel::geom_text_repel(force = force, fontface = 'bold', size = 3,
                                       segment.color = 'grey50', segment.size = 0.3)
    }
    if(!is.null(filename)){
      ggsave(plot=p, filename=filename, width = width, height = height, units = "in", ...)
    }
    return(p)
}
