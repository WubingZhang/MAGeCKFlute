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
  if(is.character(sgRNA_summary) && file.exists(sgRNA_summary)){
    dd = read.table(file = sgRNA_summary, header = TRUE, stringsAsFactors = FALSE)
  }else if(is.data.frame(sgRNA_summary) &&
           all(c("sgrna", "Gene", "LFC", "FDR") %in% colnames(sgRNA_summary))){
    dd = sgRNA_summary
  }else{
    stop("The parameter sgRNA_summary is below standard!")
  }
  dd = dd[, c("sgrna", "Gene", "LFC", "FDR")]
  rownames(dd) = dd$EntrezID
  return(dd)
}
