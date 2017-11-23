#' Gene ID conversion between ENTREZID and SYMBOL
#'
#' Gene ID conversion between ENTREZID and SYMBOL
#'
#' @docType methods
#' @name TransGeneID
#' @rdname TransGeneID
#' @aliases transGeneID
#'
#' @param genes a character vector, input genes to be converted.
#' @param fromType the input ID type, "SYMBOL" (default) or "ENTREZID".
#' @param toType the output ID type, "SYMBOL" or "ENTREZID" (default).
#' @param organism a character, specifying organism, such as "hsa" or "Human"(default), and "mmu" or "Mouse"
#'
#' @return a character vector, named by unique input gene ids.
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
#' @seealso \code{\link[pathview]{id2eg}}
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
  idx = duplicated(GeneID_Convert[,fromType])
  convert = GeneID_Convert[!idx,toType]
  names(convert) = GeneID_Convert[!idx,fromType]
  gene_after=convert[tmpGene]
  names(gene_after)=genes

  return(gene_after)
}
