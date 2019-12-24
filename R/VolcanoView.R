#' Volcano View
#'
#' Volcano plot
#'
#' @docType methods
#' @name VolcanoView
#' @rdname VolcanoView
#'
#' @param df Data frame
#' @param x Colname of df specifying x-axis in Volcanno figure, 'logFC' (default).
#' @param y Colname of df specifying y-axis in Volcanno figure, 'adj.P.Val' (default),
#' which will be plot after log10 transformation.
#' @param Label Colname of df specifying labeled terms in Volcanno figure.
#' @param top Interger, the number of top significant terms to be labeled.
#' @param topnames Character vector, indicating interested terms to be labeled.
#' @param x_cutoff Cutoff of x-axis.
#' @param y_cutoff Cutoff of y-axis.
#' @param mycolour A color vector, specifying colors of non-significant, significant up and down-regulated genes.
#' @param alpha Parameter in ggplot.
#' @param force Parameter for geom_text_repel.
#' @param main Title of volcano figure.
#' @param xlab Label of x-axis in figure.
#' @param ylab Label of y-axis in figure.
#' @param filename Figure file name to create on disk. Default filename="NULL",
#' which means don't save the figure on disk.
#' @param width Width of figure.
#' @param height Height of figure.
#' @param ... Other available parameters in ggsave.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#' @examples
#' data(rra.gene_summary)
#' rra = ReadRRA(rra.gene_summary)
#' VolcanoView(rra, x = "Score", y = "FDR", Label = "id")
#'
#' @import ggrepel
#' @export

VolcanoView <- function(df, x = "logFC", y = "adj.P.Val",
                        Label = NA, top = 5, topnames = NULL,
                        x_cutoff = log2(1.5), y_cutoff = 0.05,
                        mycolour = c("gray80", "#e41a1c", "#377eb8"),
                        alpha = 0.6, force = 0.1, main = NULL,
                        xlab = "Log2 Fold Change", ylab = "-Log10(Adjust.P)",
                        filename = NULL, width = 4, height = 2.5,
                        ...){
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
    p = ggplot(gg, aes(x=gg[,x], y=gg[,y], label=Label))
    p = p + geom_point(aes(fill=group), shape = 21, alpha=alpha, show.legend = FALSE)
    p = p + geom_point(aes(colour=color), shape = 21, alpha=alpha, show.legend = FALSE)
    p = p + scale_color_manual(values=mycolour)
    p = p + scale_fill_manual(values=mycolour)
    p = p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                  plot.title = element_text(hjust = 0.5, size=16),
                  axis.text = element_text(colour="gray10"))
    p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.border = element_blank(), panel.background = element_blank())
    p = p + geom_hline(yintercept = -log10(y_cutoff), linetype = "dotted")
    p = p + geom_vline(xintercept = c(-x_cutoff, x_cutoff), linetype = "dotted")
    p = p + labs(x=xlab, y=ylab, title=main)
    # p = p + annotate("text", color="#e41a1c", x = x_cutoff, y = max(gg[,y]), hjust = 0, vjust = 1,
    #                  label = paste("Up: ",dim(gg[gg$group=="up",])[1], sep=""))
    # p = p + annotate("text", color = "#377eb8", x = (-x_cutoff), y = max(gg[,y]), hjust = 1, vjust = 1,
    #                  label = paste("Down: ", dim(gg[gg$group=="down",])[1], sep=""))
    if(!(top==0 & is.null(topnames)))
      p = p + ggrepel::geom_text_repel(force = force, fontface = 'bold', size = 3,
                                       segment.color = 'grey50', segment.size = 0.3)
    p = p + theme(legend.position = "none")

    if(!is.null(filename)){
      ggsave(plot=p, filename=filename, width = width, height = height, units = "in", ...)
    }
    return(p)
}

