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
#' @param y Colname of df specifying y-axis in Volcanno figure, 'adj.P.Val' (default).
#' @param Label Colname of df specifying labeled terms in Volcanno figure.
#' @param top Interger, the number of top significant terms to be labeled.
#' @param topnames Character vector, indicating interested terms to be labeled.
#' @param filename Figure file name to create on disk. Default filename="NULL",
#' which means don't save the figure on disk.
#' @param x_cutoff Cutoff of x-axis.
#' @param y_cutoff Cutoff of y-axis.
#' @param main Title of volcano figure.
#' @param xlab Label of x-axis in figure.
#' @param ylab Label of y-axis in figure.
#' @param ... Other available parameters in ggsave.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#' @examples
#' data(rra.gene_summary)
#' rra = ReadRRA(rra.gene_summary)
#' VolcanoView(rra, x = "LFC", y = "FDR", Label = "Official")
#'
#' @import ggrepel
#' @export

VolcanoView <- function(df, x = "logFC", y = "adj.P.Val", Label = NA, top = 5,
                        topnames = NULL, filename = NULL, x_cutoff = log2(1.5), y_cutoff = 0.05,
                        main = NULL, xlab = "Log2 Fold Change", ylab = "-Log10(Adjust.P)", ...){
    requireNamespace("ggrepel", quietly=TRUE) || stop("need ggrepel package")
    gg = df[, c(x, y)]
    gg$group="no"
    gg$group[gg[,x]>x_cutoff & gg[,y]<y_cutoff] = "up"
    gg$group[gg[,x]< -x_cutoff & gg[,y]<y_cutoff] = "down"

    gg[, y] = -log10(gg[, y])
    if(!(top==0 & is.null(topnames))){
      gg$Label = rownames(gg)
      if(!is.na(Label)) gg$Label = df[, Label]
      gg = gg[order(abs(gg[,x]), gg[,y], decreasing = TRUE), ]
      idx1 = idx2 = c()
      if(top>0){
        idx1 = which(gg$group=="up")[1:min(top, sum(gg$group=="up"))]
        idx2 = which(gg$group=="down")[1:min(top, sum(gg$group=="down"))]
      }
      idx = unique(c(idx1, idx2, which(gg$Label %in% topnames)))
      gg$Label = as.character(gg$Label)
      gg$Label[setdiff(1:nrow(gg), idx)] = ""
      gg$Label = factor(gg$Label, levels = setdiff(unique(gg$Label), ""))
    }
    mycolour=c("no"="gray80",  "up"="#e41a1c","down"="#377eb8")
    #=========
    p = ggplot(gg, aes(x=gg[,x], y=gg[,y], colour=group, fill=group))
    p = p + geom_jitter(position = "jitter", show.legend = FALSE, alpha=0.8, size = 0.5)
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
      p = p + ggrepel::geom_text_repel(aes(x=gg[idx,x],y=gg[idx,y], label = Label), data=gg[idx,],
                                       fontface = 'bold', size = 2.5,
                                       box.padding = unit(0.4, "lines"), segment.color = 'grey50',
                                       point.padding = unit(0.3, "lines"), segment.size = 0.3)
    p = p + scale_color_manual(values=mycolour)
    p = p + scale_fill_manual(values=mycolour)
    p = p + theme(legend.position = "none")

    if(!is.null(filename)){
      ggsave(plot=p, filename=filename, units = "in", ...)
      saveRDS(df, gsub(".png|.pdf|.jpg|.tiff", ".rds", filename))
    }
    return(p)
}

