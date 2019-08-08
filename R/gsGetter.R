#' Extract pathway annotation from GMT file.
#'
#' @docType methods
#' @name gsGetter
#' @rdname gsGetter
#'
#' @param gmtpath The path to customized gmt file.
#' @param type Geneset category for testing.
#' @param limit A two-length vector (default: c(3, 50)), specifying the minimal and
#' maximal size of gene sets for enrichent analysis.
#' @param organism 'hsa' or 'mmu'.
#'
#' @return A three-column data frame.
#'
#' @author Wubing Zhang
#'
#' @examples
#' gene2path = gsGetter()
#' head(gene2path)
#'
#' @export
#'
gsGetter <- function(gmtpath = NA, type = "All", limit = c(0, Inf), organism = 'hsa'){
  if(is.na(gmtpath)){
    msigdb = file.path(system.file("extdata", package = "MAGeCKFlute"),
                       paste0(organism, "_msig_entrez.gmt.gz"))
    gene2path = ReadGMT(msigdb, limit = limit)
  }else{
    gene2path = ReadGMT(gmtpath, limit = limit)
  }
  type = unlist(strsplit(type, "\\+"))
  if(any(c("c3","c6","c7") %in% type) & (organism=='hsa')){
    msigdb = file.path(system.file("extdata", package = "MAGeCKFlute"),
                       paste0(organism, "_msig_entrez_2.gmt.gz"))
    gene2path = rbind(gene2path, ReadGMT(msigdb, limit = limit))
  }
  if("c2" %in% type) type = c(type, "KEGG", "REACTOME", "BIOCARTA", "PID")
  if("c3" %in% type) type = c(type, "GOBP", "GOCC", "GOMF")
  names(gene2path) = c("Gene","PathwayID", "PathwayName")
  gene2path$PathwayName = paste0(toupper(substr(gene2path$PathwayName, 0, 1)),
                                 tolower(substr(gene2path$PathwayName, 2,
                                                nchar(gene2path$PathwayName))))
  if(!"All" %in% type){
    idx = toupper(gsub("_.*", "", gene2path$PathwayID)) %in% toupper(type)
    if(sum(idx)>0) gene2path = gene2path[idx, ]
  }
  gene2path$Gene = as.character(gene2path$Gene)
  gene2path = gene2path[!is.na(gene2path$Gene), ]
  return(gene2path)
}
