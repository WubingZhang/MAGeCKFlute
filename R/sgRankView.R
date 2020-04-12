#' View sgRNA rank.
#'
#' @docType methods
#' @name sgRankView
#' @rdname sgRankView
#'
#' @param df A data frame, which contains columns of 'sgrna', 'Gene', and 'LFC'.
#' @param gene Character vector, specifying genes to be plotted.
#' @param top Integer, specifying number of top genes to be plotted.
#' @param bottom Integer, specifying number of bottom genes to be plotted.
#' @param neg_ctrl A vector specifying negative ctrl genes.
#' @param binwidth A numeric value specifying the bar width.
#' @param interval A numeric value specifying the interval length between each bar.
#' @param bg.col A character value specifying the background color.
#' @param filename Figure file name to create on disk. Default filename="NULL", which means no output.
#' @param width As in ggsave.
#' @param height As in ggsave.
#' @param ... Other available parameters in function 'ggsave'.
#'
#' @return An object created by \code{ggplot}.
#'
#' @author Yihan Xiao
#' @examples
#' file2 = file.path(system.file("extdata", package = "MAGeCKFlute"),
#'                   "testdata/rra.sgrna_summary.txt")
#' sgrra = ReadsgRRA(file2)
#' sgRankView(sgrra)
#'
#' @import ggplot2
#' @importFrom stats aggregate
#' @export
#'
sgRankView <- function(df, gene = NULL,
                       top = 3, bottom = 3,
                       neg_ctrl = NULL,
                       binwidth = 0.3,
                       interval = 0.1,
                       bg.col = "gray90",
                       filename = NULL,
                       width = 5, height = 3.5,
                       ...){
  if(!all(c("sgrna", "Gene", "LFC") %in% colnames(df)))
    stop("Make sure your data contains columns of 'sgrna', 'Gene', and 'LFC' ...")
  df = as.data.frame(df, stringsAsFactors = FALSE)
  df = df[order(df$LFC), ]
  df$Gene = as.character(df$Gene)
  tmp = stats::aggregate(df$LFC, by = list(df$Gene), median)
  colnames(tmp) = c("Gene", "mid")
  tmp = tmp[order(tmp$mid), ]
  if(top>0){
    idx = max((nrow(tmp)-top+1), 1)
    gene = c(gene, tmp$Gene[idx:nrow(tmp)])
  }
  if(bottom>0){
    gene = c(gene, tmp$Gene[1:min(bottom, nrow(tmp))])
  }
  gene = unique(gene)

  subdf = df[df$Gene%in%gene, ]
  if(nrow(subdf)<2) return(ggplot())
  subdf$Gene = factor(subdf$Gene, levels = gene)
  subdf = subdf[order(subdf$Gene), ]
  subdf$index = rep(1:length(gene), as.numeric(table(subdf$Gene)[gene]))
  subdf$yend <- (binwidth+interval)*subdf$index-interval
  subdf$y <- (binwidth+interval)*(subdf$index-1)
  color <- c(rep("pos",dim(subdf)[1]))
  color[which(subdf[,3]<0)] <- "neg"
  subdf$color <- color
  subdf = subdf[, c("sgrna", "Gene", "LFC", "y", "yend", "color", "index")]

  #set the scale of x-axis
  a <- -Inf
  b <- Inf

  #bgcor
  if(is.na(bg.col)){bg.col<-"white"}
  bindex <- as.vector(sapply(seq(1,max(subdf$index),1),function(x){rep(x,4)}))
  bgcol <- data.frame(as.vector(bindex))
  bgcol$color <- c(rep("bg",length(bindex)))
  colnames(bgcol)<- c("id","value")
  bgcol$x <-  rep(c(a,b,b,a),max(subdf$index))
  bgcol$y <-as.vector(sapply(seq(1,max(subdf$index),1), function(x){
    c((interval + binwidth)*(x-1), (interval + binwidth)*(x-1),
      (interval + binwidth)*x-interval,(interval + binwidth)*x-interval)
    }))

  #background
  if(!is.null(neg_ctrl)){
    neggene = rep(df[df$Gene %in% neg_ctrl, "Gene"], max(subdf$index))
    negsgrna = rep(df[df$Gene %in% neg_ctrl, "sgrna"], max(subdf$index))
    background <- data.frame(sgrna = as.vector(negsgrna), Gene = as.vector(neggene))
    background$LFC <- rep(df[df$Gene %in% neg_ctrl,3], max(subdf$index))
    seq <- as.vector(sapply(seq(1,max(subdf$index),1), function(x){
      rep(x,length(df[df$Gene %in% neg_ctrl,2]))}))
    background$y <- (binwidth+interval)*(seq-1)
    background$yend <-(binwidth+interval)*seq-interval
    background$color <-rep("tbg",length(neggene))
    background$index = 0
  }

  #depict
  cols <- c("pos"="#e41a1c","neg"="#377eb8", "tbg" = 608, "black"="black")
  p = ggplot()
  p = p + geom_polygon(aes_string("x", "y", fill="value", group="id"), color="gray20", data=bgcol)
  if(!is.null(neg_ctrl))
    p = p + geom_segment(aes_string("LFC", "y", xend = "LFC", yend = "yend", color = "color"), data = background)
  p = p + geom_segment(aes_string("LFC", "y", xend = "LFC", yend = "yend", color = "color"), data = subdf)
  p = p + scale_color_manual(values = cols)
  p = p + scale_fill_manual(values = c("bg"= bg.col))
  # p = p + scale_x_continuous(expand = c(0, 0))
  p = p + scale_y_continuous(breaks = bgcol$y[seq(1, nrow(bgcol), 4)] + binwidth/2,
                             labels = gene, expand = c(0, 0))
  p = p + labs(x = "Log2(Fold change)", y = NULL)
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
