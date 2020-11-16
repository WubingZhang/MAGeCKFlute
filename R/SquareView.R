#' Scatter plot of 9-Square
#'
#' Plot a scatter plot with Control as x-axis and Treatment as y-axis,
#' and color condition specific genes.
#'
#' @docType methods
#' @name SquareView
#' @rdname SquareView
#' @aliases squareview
#'
#' @param df Data frame, including columns of \emph{ctrlname} and \emph{treatname}, with Gene Symbol as rowname.
#' @param ctrlname A character, specifying the names of control samples.
#' @param treatname A character, specifying the name of treatment samples.
#'
#' @param label An integer or a character specifying the column used as the label, default value is 0 (row names).
#' @param label.top Boolean, whether label the top selected genes, default label the top 10 genes in each group.
#' @param top Integer, specifying the number of top selected genes to be labeled. Default is 5.
#' @param genelist Character vector, specifying labeled genes.
#'
#' @param x_cut An one or two-length numeric vector, specifying the cutoff used for x-axis.
#' @param y_cut An one or two-length numeric vector, specifying the cutoff used for y-axis.
#' @param slope A numberic value indicating slope of the diagonal cutoff.
#' @param intercept A numberic value indicating intercept of the diagonal cutoff.
#' @param auto_cut Boolean or numeric, specifying how many standard deviation will be used as cutoff.
#' @param auto_cut_x Boolean or numeric, specifying how many standard deviation will be used as cutoff on x-axis.
#' @param auto_cut_y Boolean or numeric, specifying how many standard deviation will be used as cutoff on y-axis
#' @param auto_cut_diag Boolean or numeric, specifying how many standard deviation will be used as cutoff on diagonal.
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
#' file3 = file.path(system.file("extdata", package = "MAGeCKFlute"),
#' "testdata/mle.gene_summary.txt")
#' dd = ReadBeta(file3)
#' SquareView(dd, ctrlname = "dmso", treatname = "plx", label = "Gene")
#'
#'
#' @importFrom ggrepel geom_text_repel
#' @export
#'
#'

SquareView<-function(df,
                     ctrlname = "Control",
                     treatname = "Treatment",
                     label = 0,
                     label.top = TRUE,
                     top = 5,
                     genelist = c(),
                     x_cut = NULL, y_cut = NULL,
                     slope = 1, intercept = NULL,
                     auto_cut = FALSE,
                     auto_cut_x = auto_cut,
                     auto_cut_y = auto_cut,
                     auto_cut_diag = auto_cut,
                     groups = c("midleft", "topcenter", "midright", "bottomcenter"),
                     groupnames = paste0("Group", 1:length(groups)),
                     main = NULL, filename = NULL, width = 6, height = 4, ...){
  requireNamespace("ggrepel", quietly=TRUE) || stop("need ggrepel package")

  if(!all(c(ctrlname, treatname) %in% colnames(df))){
    stop("Sample names are not consistent with column names of df.")
  }
  if(label==0) df$Gene = rownames(df)
  else df$Gene = df[, label]
  df$x = rowMeans(df[, ctrlname, drop= FALSE])
  df$y = rowMeans(df[, treatname, drop= FALSE])
  df$diff = df$y-df$x
  ## Compute the cutoff used for each dimension. ##
  if(length(x_cut)==0)
    x_cut = c(-CutoffCalling(df$x, 2), CutoffCalling(df$x, 2))
  if(length(y_cut)==0)
    y_cut = c(-CutoffCalling(df$y, 2), CutoffCalling(df$y, 2))
  if(length(intercept)==0)
    intercept = c(-CutoffCalling(df$y-slope*df$x, 2),
                  CutoffCalling(df$y-slope*df$x, 2))
  if(length(x_cut)==1)
    x_cut = sort(c(-x_cut, x_cut))
  if(length(y_cut)==1)
    y_cut = sort(c(-y_cut, y_cut))
  if(length(intercept)==1)
    intercept = sort(c(-intercept, intercept))
  ## Update the cutoff when user set the auto_cut option
  if(auto_cut_x)
    x_cut = c(-CutoffCalling(df$x, auto_cut_x),
              CutoffCalling(df$x, auto_cut_x))
  if(auto_cut_y)
    y_cut = c(-CutoffCalling(df$y, auto_cut_y),
              CutoffCalling(df$y, auto_cut_y))
  if(auto_cut_diag)
    intercept = c(-CutoffCalling(df$y-slope*df$x, auto_cut_diag),
                  CutoffCalling(df$y-slope*df$x, auto_cut_diag))

  y_min=df$x+intercept[1]; y_max=df$x+intercept[2]
  idx0 = df$y<y_max & df$y>y_min

  ## Annotate the groups ##
  df$group="Others"
  idx1 = df$x < x_cut[1]
  idx2 = df$x > x_cut[2]
  idx3 = df$y < y_cut[1]
  idx4 = df$y > y_cut[2]
  df$group[idx1&idx3] = "bottomleft"
  df$group[idx1&idx4] = "topleft"
  df$group[idx2&idx3] = "bottomright"
  df$group[idx2&idx4] = "topright"
  df$group[!idx1&!idx2&idx3] = "bottomcenter"
  df$group[!idx1&!idx2&idx4] = "topcenter"
  df$group[!idx3&!idx4&idx1] = "midleft"
  df$group[!idx3&!idx4&idx2] = "midright"
  df$group[!(df$group %in% groups)] = "Others"
  df$group[idx0] = "Others"
  mycolour=c("#377eb8", "#ff7f00", "#a65628", "#4daf4a", "#005CB7",
             "#e41a1c", "#984ea3", "#f781bf", "gray80")
  names(mycolour) = c("topleft", "topcenter", "topright", "midleft", "midright",
                      "bottomleft", "bottomcenter", "bottomright", "Others")
  names(groupnames) = groups

  ## Label top gene names ##
  df$text = df$Gene
  df$rank = top + 1
  idx = df$group=="topleft"
  df$rank[idx] = rank((df$x-df$y)[idx])
  idx = df$group=="topcenter"
  df$rank[idx] = rank(-df$y[idx])
  idx = df$group=="topright"
  df$rank[idx] = rank((-df$x-df$y)[idx])
  idx = df$group=="midleft"
  df$rank[idx] = rank((df$x)[idx])
  idx = df$group=="midright"
  df$rank[idx] = rank((-df$x)[idx])
  idx = df$group=="bottomleft"
  df$rank[idx] = rank((df$x+df$y)[idx])
  idx = df$group=="bottomcenter"
  df$rank[idx] = rank(df$y[idx])
  idx = df$group=="bottomright"
  df$rank[idx] = rank((df$y-df$x)[idx])

  df$text[df$rank>top & !(df$Gene %in% genelist)] = NA
  df$group=factor(df$group, levels = c(groups, "Others"))

  ## Limit panel region ##
  tmp = ifelse(label.top, 0.1, 0.05)
  x_min = round(min(df$x[df$group != "Others"]),2) - tmp
  x_max = round(max(df$x[df$group != "Others"]),2) + tmp
  y_min = round(min(df$y[df$group != "Others"]),2) - tmp
  y_max = round(max(df$y[df$group != "Others"]),2) + tmp
  # idx1 = (x_min<=df$x & df$x<=x_max)
  # idx2 = (y_min<=df$y  & df$y<=y_max)
  # idx = idx1&idx2
  gg = df

  ## Plot the scatter figure ##
  label_gg = gg[!is.na(gg$text), ]
  col_label = rep("#004b84", nrow(label_gg))
  col_label[label_gg$group=="Others"]="gray60"
  p = ggplot(gg, aes_string(x="x", y="y", colour="group"))
  p = p + geom_jitter(size = 1, alpha=0.8)
  p = p + scale_color_manual(values=mycolour, labels = groupnames)
  # p = p + scale_fill_discrete(guide = FALSE)
  p = p + geom_vline(xintercept = x_cut[2],linetype = "dotted")
  p = p + geom_vline(xintercept = x_cut[1],linetype = "dotted")
  p = p + geom_hline(yintercept = y_cut[2],linetype = "dotted")
  p = p + geom_hline(yintercept = y_cut[1],linetype = "dotted")
  if(!all(intercept==0)) p = p + geom_abline(slope=1, intercept=intercept, linetype = "dotted")
  p = p + labs(x="Control", y = "Treatment")
  # p = p + guides(col = guide_legend(ncol = 3, byrow = TRUE))
  if(label.top)
    p = p + ggrepel::geom_text_repel(aes_string(x="x",y="y",label="Gene"),
                                     color=col_label, data=label_gg)
  if("topleft" %in% groups)
    p = p + annotate("text", color="red", x = x_cut[1], y=max(gg$y),
                     vjust=0, hjust = 3, label=sum(gg$group=="topleft"))
  if("topcenter" %in% groups)
    p = p + annotate("text", color="red", x = 0, y=max(gg$y),
                     vjust=0, hjust = 0, label=sum(df$group=="topcenter"))
  if("topright" %in% groups)
    p = p + annotate("text", color="red", x = x_cut[2], y=max(gg$y),
                     vjust=0, hjust = -3, label=sum(df$group=="topright"))
  if("midleft" %in% groups)
    p = p + annotate("text", color="red", x = min(gg$x), y=0,
                     vjust=0, hjust = -1, label=sum(df$group=="midleft"))
  if("midright" %in% groups)
    p = p + annotate("text", color="red", x = max(gg$x), y=0,
                     vjust=0, hjust = 1, label=sum(df$group=="midright"))
  if("bottomleft" %in% groups)
    p = p + annotate("text", color="red", x = x_cut[1], y=min(gg$y),
                     vjust=0, hjust = 3, label=sum(df$group=="bottomleft"))
  if("bottomcenter" %in% groups)
    p = p + annotate("text", color="red", x = 0, y=min(gg$y),
                     vjust=0, hjust = 0, label=sum(df$group=="bottomcenter"))
  if("bottomright" %in% groups)
    p = p + annotate("text", color="red", x = x_cut[2], y=min(gg$y),
                     vjust=0, hjust = -3, label=sum(df$group=="bottomright"))
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
  # p = p + xlim(x_min, x_max) + ylim(y_min, y_max)

  if(!is.null(filename)){
      write.table(p$data, gsub("\\....$", ".txt", filename),
                sep = "\t", quote = FALSE, row.names = FALSE)
      ggsave(filename, p, units="in", width=width, height=height, ...)
  }
  return(p)
}
