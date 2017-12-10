#' Normalize gene beta scores
#'
#' Two normalization methods are available. \code{cell_cycle} method normalizes gene beta scores
#'  based on positive control genes in CRISPR screening. \code{loess} method normalizes gene
#'  beta scores using loess.
#'
#' @docType methods
#' @name NormalizeBeta
#' @rdname NormalizeBeta
#' @aliases normalizebeta
#'
#' @param beta data frame, which has columns of 'Gene', and \code{samples} beginning from the third columns.
#' @param samples character vector, specifying the samples in \code{beta} to be normalized.
#' @param method character, one of 'cell_cycle'(default) and 'loess'.
#' @param posControl a file path or a character vector, specifying positive control genes used
#' for cell cycle normalization.
#' @param minus numeric, scale for cell cycle normalization. Between 0 and 1.
#'
#' @return A data frame with same format as input data \code{beta}.
#'
#' @details In CRISPR screens, cells treated with different conditions (e.g., with or without
#' drug) may have different proliferation rates. So we defined a list of core essential genes,
#' which is equally negatively selected between samples with different proliferation rate.
#' Normalization of gene beta scores is performed using these essential genes. \code{cell_cycle}
#' in MAGeCKFlute normalizes the beta scores of all genes based on the median beta score of essential genes.
#' After normalization, the beta scores are comparable across samples. \code{loess} is another
#' optional normalization method, which is used to normalize array data before.
#'
#' @author Wubing Zhang
#'
#' @note See the vignette for an example of NormalizeBeta.
#' Note that the source code of \code{NormalizeBeta} is very simple.
#' The source can be found by typing \code{MAGeCKFlute:::NormalizeBeta}
#' or \code{getMethod("NormalizeBeta")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/NormalizeBeta.R}
#' Users should find it easy to customize this function.
#'
#'
#' @examples
#' data(MLE_Data)
#' # Read beta score from gene summary table in MAGeCK MLE results
#' dd = ReadBeta(MLE_Data, organism="hsa")
#' #Cell Cycle normalization
#' dd_essential = NormalizeBeta(dd, method="cell_cycle")
#' head(dd_essential)
#'
#' #Optional loess normalization
#' dd_loess = NormalizeBeta(dd, method="loess")
#' head(dd_loess)
#'
#' @importFrom Biobase rowMedians
#' @importFrom affy normalize.loess
#'
#' @export

#===normalize function=====================================
NormalizeBeta <- function(beta, samples=NULL,
                          method="cell_cycle", posControl=NULL, minus=0.6){
  loginfo("Normalize beta scores ...")
  requireNamespace("Biobase", quietly=TRUE) || stop("need Biobase package")
  if(is.null(samples)) samples = colnames(beta)[3:ncol(beta)]

  if(method=="cell_cycle"){
    if(!is.null(posControl) && class(posControl)=="character" && file.exists(posControl)[1]){
      tmp = read.table(posControl, sep = "\t", header = FALSE)
      posControl = as.character(unlist(tmp))
    }else{
      data(Core_Essential)
      posControl=Core_Essential
    }
    idx = which(beta$Gene %in% posControl)
    normalized = as.matrix(beta[,samples])
    mid = rowMedians(t(normalized[idx,]))
    mid = abs(mid - minus)
    normalized = t(t(normalized) / mid)
  }
  if(method=="loess"){
    requireNamespace("affy", quietly=TRUE) || stop("need affy package")
    normalized = as.matrix(beta[,samples])
    normalized = normalize.loess(normalized,log.it = FALSE, verbose=FALSE)
  }
  beta[,samples] = normalized

  return(beta)
}
