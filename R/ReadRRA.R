#' Read MAGeCK-RRA data
#'
#' Read pvalue of gene selection from file or data frame
#'
#' @docType methods
#' @name ReadRRA
#' @rdname ReadRRA
#' @aliases readrra
#'
#' @param gene_summary a file path or a data frame, which has three columns named 'id', 'neg.fdr' and 'pos.fdr'.
#' @param organism character, KEGG species code, or the common species name, used to determine
#' the gene annotation package. For all potential values check: data(bods); bods. Default org="hsa",
#' and can also be "human" (case insensitive).
#'
#' @return a data frame including four columns, named "Official", "neg.fdr", "pos.fdr" and "ENTREZID".
#'
#' @author Wubing Zhang
#'
#' @note
#' The source can be found by typing \code{MAGeCKFlute:::ReadRRA}
#' or \code{getMethod("ReadRRA")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/ReadRRA.R}
#' Users should find it easy to customize this function.
#'
#'
#' @examples
#' data(RRA_Data)
#' dd.rra = ReadRRA(RRA_Data, organism="hsa")
#' head(dd.rra)
#'
#' @import pathview
#'
#' @export
#'


#===read gene summary file=====================================
ReadRRA <- function(gene_summary, organism="hsa"){
  loginfo("Read gene summary file ...")
  if(class(gene_summary)=="character" && file.exists(gene_summary)){
    dd=read.table(file=gene_summary,header= TRUE)
  }else if(class(gene_summary)=="data.frame" &&
           all(c("id", "neg.fdr","pos.fdr")%in%colnames(gene_summary))){
    dd=gene_summary
  }else{
    stop("The parameter gene_summary is below standard!")
  }
  idx=grepl("Zhang_", dd$id, ignore.case = TRUE)
  dd=dd[!idx,]
  idx=grepl("^CTRL", dd$id, ignore.case = TRUE)
  dd=dd[!idx,]

  dd=dd[,c("id","neg.fdr","pos.fdr")]
  dd$ENTREZID = suppressMessages(id2eg(dd$id, "SYMBOL", org = organism)[, "ENTREZID"])
  idx=is.na(dd$ENTREZID)
  dd=dd[!idx,]
  colnames(dd) = c("Official", "neg.fdr", "pos.fdr", "ENTREZID")

  return(dd)
}
