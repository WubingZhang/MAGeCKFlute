#' Determine the gene annotation package.
#'
#' Determine the gene annotation package. for specific organism
#'
#' @docType methods
#' @name getOrg
#' @rdname getOrg
#'
#' @param organism Character, KEGG species code, or the common species name, used to determine
#' the gene annotation package. For all potential values check: data(bods); bods. Default org="hsa",
#' and can also be "human" (case insensitive).
#'
#' @return A list containing three elements:
#' \item{organism}{species}
#' \code{pkg}{annotation package name}
#' \code{Symbol_Entrez}{a data frame, mapping between gene symbol and entrez id}
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
  requireNamespace("clusterProfiler", quietly=TRUE) || stop("need clusterProfiler package")
  data(bods)
  res=list()
  ##======
  # Get the mapping from organism to package
  ridx=grep(tolower(organism), tolower(bods[,2:3])) %% nrow(bods)
  if(length(ridx)==0) stop("Wrong organism name!")
  res$org = bods[ridx,3]
  res$pkg = bods[ridx,1]

  # Get the mapping from organism to alias2eg
  tmp = strsplit(res$pkg, "\\.")[[1]]
  tmp[3] = "egALIAS2EG"
  tmp = tmp[-4]
  alias = paste(tmp, collapse = ".")
  #==============
  # Get mapping
  require(res$pkg, character.only = TRUE)
  tmpFile = file.path(temp_dir(), paste0("map_symbol_entrez_", res$org))
  if(!file.exists(tmpFile)){
    aliasMapping <- as.data.frame(toTable(get(alias)))
    tmp=suppressWarnings(suppressMessages(as.data.frame(
      bitr(aliasMapping$gene_id, fromType="ENTREZID", toType="SYMBOL", OrgDb=res$pkg))))
    aliasMapping = merge(aliasMapping,tmp,by.x="gene_id",by.y="ENTREZID")
    aliasMapping = data.frame(gene_id=aliasMapping$gene_id,
                              SYMBOL=c(aliasMapping$SYMBOL, aliasMapping$alias_symbol))
    tmp=suppressWarnings(suppressMessages(as.data.frame(
      bitr(aliasMapping$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb=res$pkg))))
    tmp1 = merge(aliasMapping,tmp,by.x="SYMBOL",by.y="SYMBOL")
    mapping = tmp1[,c("SYMBOL","ENTREZID")]
    idx1 = duplicated(mapping$SYMBOL)
    idx2 = duplicated(mapping$ENTREZID)
    mapping = mapping[!(idx1&idx2),]
    mapping$SYMBOL = toupper(mapping$SYMBOL)
    mapping$ENTREZID = as.character(mapping$ENTREZID)
    res$Symbol_Entrez = mapping
    write.table(mapping, tmpFile, sep = "\t", row.names = FALSE)
  }else{
    res$Symbol_Entrez = read.table(tmpFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  }
  return(res)
}
