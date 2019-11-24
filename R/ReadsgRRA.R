#' Read sgRNA summary in MAGeCK-RRA results
#'
#' @docType methods
#' @name ReadsgRRA
#' @rdname ReadsgRRA
#'
#' @param sgRNA_summary A file path or a data frame of sgRNA summary data.
#'
#' @return A data frame.
#'
#' @author Wubing Zhang
#'
#' @examples
#' data(rra.sgrna_summary)
#' sgrra = ReadsgRRA(rra.sgrna_summary)
#' head(sgrra)
#'
#' @export
#'
ReadsgRRA <- function(sgRNA_summary){
  message(Sys.time(), " # Read sgRNA summary file ...")
  if(is.null(dim(sgRNA_summary))){
    sgRNA_summary = read.table(file = sgRNA_summary, sep = "\t", header = TRUE, quote = "",
                    comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
  }
  dd = sgRNA_summary[, c("sgrna", "Gene", "LFC", "FDR")]
  return(dd)
}
