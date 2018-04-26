#' Gene ID conversion between ENTREZID and SYMBOL
#'
#' Gene ID conversion between ENTREZID and SYMBOL
#'
#' @docType methods
#' @name TransGeneID
#' @rdname TransGeneID
#' @aliases transGeneID
#'
#' @param genes A character vector, input genes to be converted.
#' @param fromType The input ID type, one of "SYMBOL" (default), "ENTREZID",
#' "ENSEMBL", "HGNC", "FULLNAME" and "TYPE".
#' @param toType The output ID type, one of "SYMBOL", "ENTREZID" (default),
#' "ENSEMBL", "HGNC", "FULLNAME" and "TYPE".
#' @param organism A character, specifying organism, one of "hsa"(or 'Human'),
#' "mmu"(or 'Mouse'), "bta", "cfa", "ptr", "rno", and "ssc".
#'
#' @return A character vector, named by unique input gene ids.
#'
#' @author Wubing Zhang
#'
#' @note
#' The source can be found by typing \code{MAGeCKFlute:::TransGeneID}
#' or \code{getMethod("TransGeneID")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/TransGeneID.R}
#' Users should find it easy to customize this function.
#'
#' @seealso \code{\link[pathview]{eg2id}}
#'
#' @examples
#' data(MLE_Data)
#' TransGeneID(MLE_Data$Gene[1:10], fromType="SYMBOL", toType="ENTREZID", organism="hsa")
#'
#' @export


TransGeneID <- function(genes, fromType="SYMBOL", toType="ENTREZID", organism="hsa"){
  if(is.null(genes) | length(genes)==0){
    return(c())
  }
  genes = as.character(genes)
  tmpGene = toupper(genes)
  #=====================
  GeneID_Convert = getOrg(organism)$Symbol_Entrez
  idx = duplicated(GeneID_Convert[, fromType])
  convert = GeneID_Convert[!idx, toType]
  names(convert) = GeneID_Convert[!idx,fromType]
  gene_after=convert[tmpGene]
  names(gene_after)=genes

  return(gene_after)
}
