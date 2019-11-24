#' Read gene beta scores
#'
#' Read gene beta scores from file or data frame
#'
#' @docType methods
#' @name ReadBeta
#' @rdname ReadBeta
#' @aliases readbeta
#'
#' @param gene_summary A file path or a data frame, which has columns of 'Gene' and beta score of samples.
#'
#' @return A data frame, in which the first column is ENTREZID, and the later columns are beta score for each samples.
#'
#' @author Wubing Zhang
#'
#' @examples
#' data(mle.gene_summary)
#' dd = ReadBeta(mle.gene_summary)
#' head(dd)
#'
#' @export
ReadBeta <- function(gene_summary){
  message(Sys.time(), " # Read gene summary file ...")
  if(is.null(dim(gene_summary))){
    gene_summary = read.table(file = gene_summary, sep = "\t", header = TRUE, quote = "",
                              comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
  }
  dd = gene_summary[, c(1,seq(3,ncol(gene_summary),6))]
  names(dd) = gsub(".beta", "", names(dd))
  return(dd)
}
