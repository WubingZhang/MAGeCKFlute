#' Rank plot
#'
#' Draw the score and rank of genes on a scatter plot.
#'
#' @docType methods
#' @name RankView
#' @rdname RankView
#' @aliases rankview
#'
#' @param rankdata A numeric vector, with gene as names.
#' @param genelist A character vector, specifying genes to be labeled.
#' @param decreasing Boolean, specifying the order of genes to plot.
#' @param top Integer, specifying number of positive genes to be labeled.
#' @param bottom Integer, specifying number of negative genes to be labeled.
#' @param cutoff One numeric value indicating the fold of standard deviation used as cutoff;
#' two number vector, such as c(-1, 1), specifying the exact cutoff for selecting top genes.
#' @param main A character, specifying title.
#' @param filename A character, specifying a file name to create on disk.
#' Set filename to be "NULL", if don't want to save the figure.
#' @param width Numeric, specifying width of figure.
#' @param height Numeric, specifying height of figure.
#' @param ... Other available parameters in the function 'geom_text_repel'.
#'
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#'
#' @examples
#' file1 = file.path(system.file("extdata", package = "MAGeCKFlute"),
#' "testdata/rra.gene_summary.txt")
#' gdata = ReadRRA(file1)
#' rankdata = gdata$Score
#' names(rankdata) = gdata$id
#' RankView(rankdata)
#'
#' @importFrom ggrepel geom_text_repel
#'
#' @export
#'

RankView <- function(rankdata, genelist = NULL, decreasing = TRUE,
                     top = 5, bottom = 5, cutoff = 2,
                     main = NULL, filename = NULL, width = 5, height = 4, ...){
  requireNamespace("ggrepel", quietly=TRUE) || stop("need ggrepel package")
  if(length(cutoff)<2){
    cutoff = CutoffCalling(rankdata, cutoff)
    cutoff = c(-cutoff, cutoff)
  }else cutoff = sort(cutoff)
  data = data.frame(Gene = names(rankdata), diff = rankdata, stringsAsFactors=FALSE)
  if(decreasing) data$Rank = rank(-data$diff) else data$Rank = rank(data$diff)
  data$group = "no"
  data$group[data$diff>cutoff[2]] = "up"
  data$group[data$diff<cutoff[1]] = "down"

  idx=(data$Rank<=bottom) | (data$Rank>(max(data$Rank)-top)) | (data$Gene %in% genelist)
  mycolour = c("no"="gray80",  "up"="#e41a1c","down"="#377eb8")
  p = ggplot(data, aes_string(x="Rank",y="diff",color="group"))
  p = p + geom_point(size = 0.5)
  # p = p + geom_hline(yintercept = cutoff, linetype = "dotted")
  if(sum(idx)>0)
    p = p + ggrepel::geom_text_repel(aes_string(label = "Gene"), data=data[idx,], size = 2.5, ...)
  p = p + scale_color_manual(values=mycolour)
  p = p + scale_fill_manual(values=mycolour)
  p = p + labs(x = "Rank", y = "Score", title=main)
  p = p + theme_bw(base_size = 14)
  p = p + theme(plot.title = element_text(hjust = 0.5))
  p = p + theme(legend.position="none")#+ylim(-1000,7000)

  if(!is.null(filename)){
    ggsave(plot=p, filename=filename, units = "in", width=width, height=height)
  }
  return(p)
}
