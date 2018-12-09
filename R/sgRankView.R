#' View the rank of sgRNA points in horizontal bars.
#'
#' @docType methods
#' @name sgRankView
#' @rdname sgRankView
#'
#' @param df A data frame, which contains columns of 'sgrna', 'Gene', and 'LFC'.
#' @param genelist Character vector, specifying genes to be plotted.
#' @param top Integer, specifying number of top genes to be plotted.
#' @param bottom Integer, specifying number of bottom genes to be plotted.
#' @param neg_ctrl A vector of negative ctrl genes.
#' @param binwidth A numeric value specifying the bar width,
#' @param line.size A numeric value specifying the size of segment lines.
#' @param filename Figure file name to create on disk. Default filename="NULL", which means no output.
#' @param width As in ggsave.
#' @param height As in ggsave.
#' @param ... Other available parameters in function 'ggsave'.
#'
#' @return An object created by \code{ggplot}.
#'
#' @author Wubing Zhang
#' @examples
#' data(rra.sgrna_summary)
#' sgrra = ReadsgRRA(rra.sgrna_summary)
#' sgRankView(sgrra)
#'
#' @import ggplot2
#' @export
#'
sgRankView <- function(df, genelist = NULL, top = 3, bottom = 3, neg_ctrl = NULL,
                       binwidth = 5, line.size = 1,
                       filename = NULL, width = 5, height = 3.5, ...){
  df = df[order(df$LFC), ]
  df$sgrna = as.character(df$sgrna)
  df$Gene = as.character(df$Gene)
  if(top>0){
    top_genes = sort(table(df$Gene[(nrow(df)-100*top+1):nrow(df)]), decreasing = TRUE)
    genelist = c(genelist, names(top_genes[1:top]))
  }
  if(bottom>0){
    bott_genes = sort(table(df$Gene[1:(bottom*100)]), decreasing = TRUE)
    genelist = c(genelist, names(bott_genes[1:bottom]))
  }
  genelist = unique(genelist)
  idx1 = which(df$Gene %in% genelist)
  gg = df[idx1, ]
  gg$group = "pos"
  neg_gene = sapply(unique(gg$Gene), function(x) median(gg$LFC[gg$Gene==x]))
  gg$group[gg$Gene%in%names(neg_gene[neg_gene<0])] = "neg"
  gg$Gene = factor(gg$Gene, levels = names(sort(neg_gene, decreasing = TRUE)))
  idx2 = df$Gene %in% neg_ctrl
  if(sum(idx2)>0){
    for(i in 1:length(unique(gg$Gene))){
      tmp_gg = df[idx2, ]
      tmp_gg$Gene = unique(gg$Gene)[i]
      tmp_gg$group = "no"
      gg = rbind.data.frame(gg, tmp_gg)
    }
  }
  gg$x_end = gg$LFC - line.size*0.02
  gg$x_end[gg$LFC<0] = gg$LFC[gg$LFC<0] + line.size*0.02

  # Plot segments
  p = ggplot(gg, aes(x = LFC, y = Gene))
  p = p + geom_segment(aes(x=min(df$LFC), xend=max(df$LFC), yend=Gene),
                       size=binwidth+5+line.size*0.5, colour="black", lineend = "butt")
  p = p + geom_segment(aes(x=min(df$LFC)+line.size*0.01,
                           xend=max(df$LFC)-line.size*0.01, yend=Gene),
                       size=binwidth+5, colour="gray80", lineend = "butt")
  p = p + geom_segment(aes(xend = x_end, yend = Gene, color = group), size=binwidth+5)
  p = p + scale_color_manual(values = c("pos"="#e41a1c","neg"="#377eb8", "no"="gray50"))
  p = p + labs(x = NULL, y = NULL)
  p = p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_blank(),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank())
  p = p + theme(legend.position = "none")

  if(!is.null(filename)){
    ggsave(plot=p, filename=filename, units = "in", width=width, height=height, ...)
  }
  return(p)
}
