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
  if(is.character(gene_summary) && file.exists(gene_summary)){
    dd=read.table(file = gene_summary, header = TRUE,
                  check.names = FALSE, stringsAsFactors = FALSE)
  }else if(is.data.frame(gene_summary) &&
          ("Gene" %in% names(gene_summary)) &&
          any(grepl(".beta", names(gene_summary)))){
    dd = gene_summary
  }else{
    stop("gene_summary is invalid!")
  }

  idx = grepl(".beta",names(dd))
  idx[1] = TRUE; dd = dd[,idx]
  names(dd) = gsub(".beta", "", names(dd))

  return(dd)
}
