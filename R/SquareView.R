#' Scatter plot of 9-Square
#'
#' Plot a scatter plot with Control beta score as x-axis and Treatment beta score as y-axis,
#' and colored treatment related genes.
#'
#' @docType methods
#' @name SquareView
#' @rdname SquareView
#' @aliases squareview
#'
#' @param beta Data frame, including columns of \emph{ctrlname} and \emph{treatname}, with Gene Symbol as rowname.
#' @param ctrlname A character, specifying the names of control samples.
#' @param treatname A character, specifying the name of treatment samples.
#'
#' @param label An integer or a character specifying the column used as the label, default value is 0 (row names).
#' @param label.top Boolean, whether label the top selected genes, default label the top 10 genes in each group.
#' @param top Integer, specifying the number of top selected genes to be labeled. Default is 5.
#' @param genelist Character vector, specifying labeled genes.
#'
#' @param x_cutoff An one or two-length numeric vector, specifying the cutoff used for x-axis.
#' @param y_cutoff An one or two-length numeric vector, specifying the cutoff used for y-axis.
#' @param intercept An one or two-length numeric vector, specifying the intercept of diagonal.
#'
#' @param groups A character vector, specifying which group to be colored. Optional groups include "topleft",
#' "topcenter", "topright", "midleft", "midright", "bottomleft", "bottomcenter", "bottomright".
#' @param groupnames A character vector, specifying group names.
#'
#' @param main As in 'plot'.
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
#' @seealso \code{\link{ScatterView}}
#'
#' @examples
#' data(mle.gene_summary)
#' # Read beta score from gene summary table in MAGeCK MLE results
#' dd = ReadBeta(mle.gene_summary, organism="hsa")
#' SquareView(dd, ctrlname = "dmso", treatname = "plx", label = "Gene")
#'
#'
#' @importFrom ggExtra ggMarginal
#' @importFrom ggrepel geom_text_repel
#' @export
#'
#'

SquareView<-function(beta, ctrlname = "Control", treatname = "Treatment",
                     label = 0, label.top = TRUE, top = 5, genelist = c(),
                     x_cutoff = c(-1,1), y_cutoff = c(-1,1), intercept = NULL,
                     groups = c("midleft", "topcenter", "midright", "bottomcenter"),
                     groupnames = paste0("Group", 1:length(groups)),
                     main = NULL, filename = NULL, width = 5, height = 4, ...){
  requireNamespace("ggExtra", quietly=TRUE) || stop("need ggExtra package")
  requireNamespace("ggrepel", quietly=TRUE) || stop("need ggrepel package")

  message(Sys.time(), " # Square plot for ", main, " beta scores ...")
  if(!all(c(ctrlname, treatname) %in% colnames(beta))){
    stop("Sample names are not consistent with column names of beta.")
  }
  if(label==0) beta$Gene = rownames(beta)
  else beta$Gene = beta[, label]
  beta$x = rowMeans(beta[, ctrlname, drop= FALSE])
  beta$y = rowMeans(beta[, treatname, drop= FALSE])
  beta$diff = beta$y-beta$x
  ## Compute the cutoff used for each dimension. ##
  if(length(x_cutoff)==1) x_cutoff = sort(c(-x_cutoff, x_cutoff))
  if(length(y_cutoff)==1) y_cutoff = sort(c(-y_cutoff, y_cutoff))
  if(length(intercept)==0) intercept = CutoffCalling(beta$diff)
  if(length(intercept)==1) intercept = sort(c(-intercept, intercept))
  y_min=beta$x+intercept[1]; y_max=beta$x+intercept[2]
  idx0 = beta$y<y_max & beta$y>y_min

  ## Annotate the groups ##
  beta$group="Others"
  idx1 = beta$x < x_cutoff[1]
  idx2 = beta$x > x_cutoff[2]
  idx3 = beta$y < y_cutoff[1]
  idx4 = beta$y > y_cutoff[2]
  beta$group[idx1&idx3] = "bottomleft"
  beta$group[idx1&idx4] = "topleft"
  beta$group[idx2&idx3] = "bottomright"
  beta$group[idx2&idx4] = "topright"
  beta$group[!idx1&!idx2&idx3] = "bottomcenter"
  beta$group[!idx1&!idx2&idx4] = "topcenter"
  beta$group[!idx3&!idx4&idx1] = "midleft"
  beta$group[!idx3&!idx4&idx2] = "midright"
  beta$group[!(beta$group %in% groups)] = "Others"
  beta$group[idx0] = "Others"
  mycolour=c("#377eb8", "#ff7f00", "#a65628", "#4daf4a", "#005CB7",
             "#e41a1c", "#984ea3", "#f781bf", "gray80")
  names(mycolour) = c("topleft", "topcenter", "topright", "midleft", "midright",
                      "bottomleft", "bottomcenter", "bottomright", "Others")
  names(groupnames) = groups

  ## Label top gene names ##
  beta$text = beta$Gene
  beta$rank = top + 1
  idx = beta$group=="topleft"
  beta$rank[idx] = rank((beta$x-beta$y)[idx])
  idx = beta$group=="topcenter"
  beta$rank[idx] = rank(-beta$y[idx])
  idx = beta$group=="topright"
  beta$rank[idx] = rank((-beta$x-beta$y)[idx])
  idx = beta$group=="midleft"
  beta$rank[idx] = rank((beta$x)[idx])
  idx = beta$group=="midright"
  beta$rank[idx] = rank((-beta$x)[idx])
  idx = beta$group=="bottomleft"
  beta$rank[idx] = rank((beta$x+beta$y)[idx])
  idx = beta$group=="bottomcenter"
  beta$rank[idx] = rank(beta$y[idx])
  idx = beta$group=="bottomright"
  beta$rank[idx] = rank((beta$y-beta$x)[idx])

  beta$text[beta$rank>top & !(beta$Gene %in% genelist)] = NA
  beta$group=factor(beta$group, levels = c(groups, "Others"))

  ## Limit panel region ##
  tmp = ifelse(label.top, 0.1, 0.05)
  x_min = round(min(beta$x[beta$group != "Others"]),2) - tmp
  x_max = round(max(beta$x[beta$group != "Others"]),2) + tmp
  y_min = round(min(beta$y[beta$group != "Others"]),2) - tmp
  y_max = round(max(beta$y[beta$group != "Others"]),2) + tmp
  idx1 = (x_min<=beta$x & beta$x<=x_max)
  idx2 = (y_min<=beta$y  & beta$y<=y_max)
  idx = idx1&idx2
  gg = beta[idx, ]

  ## Plot the scatter figure ##
  label_gg = gg[!is.na(gg$text), ]
  col_label = rep("#004b84", nrow(label_gg))
  col_label[label_gg$group=="Others"]="gray60"
  p = ggplot(gg, aes(x=x, y=y, colour=group))
  p = p + geom_jitter(size = 1, alpha=0.8)
  p = p + scale_color_manual(values=mycolour, labels = groupnames)
  # p = p + scale_fill_discrete(guide = FALSE)
  p = p + geom_vline(xintercept = x_cutoff[2],linetype = "dotted")
  p = p + geom_vline(xintercept = x_cutoff[1],linetype = "dotted")
  p = p + geom_hline(yintercept = y_cutoff[2],linetype = "dotted")
  p = p + geom_hline(yintercept = y_cutoff[1],linetype = "dotted")
  p = p + geom_abline(slope=1, intercept=intercept, linetype = "dotted")
  p = p + labs(x="Control beta score", y = "Treatment beta score")
  # p = p + guides(col = guide_legend(ncol = 3, byrow = TRUE))
  if(label.top)
    p = p + ggrepel::geom_text_repel(aes(x=x,y=y,label=Gene), color=col_label, data=label_gg)
  if("topleft" %in% groups)
    p = p + annotate("text", color="red", x = x_cutoff[1], y=max(gg$y),
                     vjust=0, hjust = 3, label=sum(gg$group=="topleft"))
  if("topcenter" %in% groups)
    p = p + annotate("text", color="red", x = 0, y=max(gg$y),
                     vjust=0, hjust = 0, label=sum(beta$group=="topcenter"))
  if("topright" %in% groups)
    p = p + annotate("text", color="red", x = x_cutoff[2], y=max(gg$y),
                     vjust=0, hjust = -3, label=sum(beta$group=="topright"))
  if("midleft" %in% groups)
    p = p + annotate("text", color="red", x = min(gg$x), y=0,
                     vjust=0, hjust = -1, label=sum(beta$group=="midleft"))
  if("midright" %in% groups)
    p = p + annotate("text", color="red", x = max(gg$x), y=0,
                     vjust=0, hjust = 1, label=sum(beta$group=="midright"))
  if("bottomleft" %in% groups)
    p = p + annotate("text", color="red", x = x_cutoff[1], y=min(gg$y),
                     vjust=0, hjust = 3, label=sum(beta$group=="bottomleft"))
  if("bottomcenter" %in% groups)
    p = p + annotate("text", color="red", x = 0, y=min(gg$y),
                     vjust=0, hjust = 0, label=sum(beta$group=="bottomcenter"))
  if("bottomright" %in% groups)
    p = p + annotate("text", color="red", x = x_cutoff[2], y=min(gg$y),
                     vjust=0, hjust = -3, label=sum(beta$group=="bottomright"))
  p = p + theme(legend.key = element_rect(fill = "transparent"))
  p = p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank())
  # p = p + theme(legend.position="none")
  p = p + theme(legend.title=element_blank())
  # p = suppressWarnings(ggExtra::ggMarginal(p, type="histogram", bins=50, fill = "gray80"))
  p$data = beta[, c("Gene", "EntrezID", "x", "y", "diff", "group", "text")]
  p$data = p$data[order(p$data$group), ]

  if(!is.null(filename)){
      write.table(p$data, gsub("\\....$", ".txt", filename),
                sep = "\t", quote = FALSE, row.names = FALSE)
      ggsave(filename, p, units="in", width=width, height=height, ...)
  }
  return(p)
}
