#' Determine the gene annotation package.
#'
#' Determine the gene annotation package. for specific organism
#'
#' @docType methods
#' @name getOrg
#' @rdname getOrg
#'
#' @param organism character, KEGG species code, or the common species name, used to determine
#' the gene annotation package. For all potential values check: data(bods); bods. Default org="hsa",
#' and can also be "human" (case insensitive).
#'
#' @return a list containing two items:
#' \item{organism}{species}
#' \code{pkg}{annotation package name}
#'
#' @author Wubing Zhang
#'
#' @note
#' The source can be found by typing \code{MAGeCKFlute:::getOrg}
#' or \code{getMethod("getOrg")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/getOrg.R}
#' Users should find it easy to customize this function.
#'
#' @examples
#' ann_pkg = getOrg("human")$pkg
#' print(ann_pkg)
#'
#' @importFrom reshape melt

#'
#' @export

getOrg <- function(organism){
  requireNamespace("reshape", quietly=TRUE) || stop("need reshape package")
  data(bods)
  res=list()
  ##======
  # Get the mapping from organism to package
  orgMap = as.data.frame(bods[,1:3], stringsAsFactors=FALSE)
  orgMap$species = tolower(orgMap$species)
  orgMap = melt(orgMap, id="package")
  if(tolower(organism) %in% orgMap$value){
    res$pkg = orgMap[match(tolower(organism), orgMap[, 3]), 1]
    res$org = as.character(orgMap[19+match(res$pkg, orgMap[20:38, 1]), 3])
  }else{
    stop("No organism found !")
  }

  return(res)

}
