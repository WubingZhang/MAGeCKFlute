#' Select signatures from candidate list (according to the consistence in most samples).
#'
#' @docType methods
#' @name Selector
#' @rdname Selector
#'
#' @param mat Data matrix, each row is candidates (genes), each column is samples.
#' @param cutoff Cutoff to define the signatures.
#' @param type Direction to select signatures.
#' @param select Proportion of samples in which signature is selected.
#'
#' @return An list containing two elements, first is selected signature and second is a ggplot object.
#'
#' @examples
#' mat = matrix(rnorm(1000*30), 1000, 30)
#' rownames(mat) = paste0("Gene", 1:1000)
#' colnames(mat) = paste0("Sample", 1:30)
#' hits = Selector(mat, select = 0.68)
#' print(hits$p)
#'
#' @import ggsci
#' @export

Selector <- function(mat, cutoff = 0, type = "<", select = 0.8){
  Num = 1:ncol(mat)
  gene_count = apply(mat, 1, function(x)
    ifelse(type==">", length(which(x>=cutoff)), length(which(x<=cutoff))))
  gg = data.frame(Num=Num, GeneNum = sapply(Num, function(x)
                            length(which(gene_count>=x))))

  text <- data.frame(x=c(gg[1,1], round(select*ncol(mat)), ncol(mat)),
                     y=c(gg[1,2], gg[gg$Num==round(select*ncol(mat)), 2], gg[ncol(mat), 2]))
  genes = names(gene_count)[gene_count>=round(select*ncol(mat))]
  #============
  p = ggplot(gg)
  p = p + geom_line(aes(x=Num,y=GeneNum),colour = "#2972A3")
  p = p + geom_vline(xintercept = round(select*ncol(mat)),linetype=2,size=0.5)
  p = p + geom_hline(yintercept = gg[gg$Num==round(select*ncol(mat)), 2],linetype=2,size=0.5)
  p = p + geom_point(aes(x=x,y=y),data=text, shape=17, size=2, color="#d66648")
  p = p + geom_text(aes(x=x,y=y),data=text, label=paste("(",text$x,", ",text$y,")",sep=""),
                    colour = "black",size=3.5)
  p = p + ylab("Gene number")+xlab("Sample number")+scale_colour_npg()
  p = p + theme(text = element_text(colour="black",size = 14),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank())

  return(list(sig=genes, p = p))
}
