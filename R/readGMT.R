#' ReadGMT
#'
#' Parse gmt file to a data.frame
#'
#' @docType methods
#' @name ReadGMT
#' @rdname ReadGMT
#'
#' @param gmtpath The path to gmt file.
#' @param limit A integer vector of length two, specifying the limit of geneset size.
#'
#' @return An data.frame, in which the first column is gene, and the second column is pathway name.
#'
#' @author Wubing Zhang
#'
#' @export
#'
ReadGMT <- function(gmtpath, limit = c(0, Inf)){
  c2 = readLines(gmtpath)
  c2_list = strsplit(c2, "\t")
  c2_pathway = t(matrix(unlist(lapply(c2_list, function(x) x[1:2])), nrow = 2))
  rownames(c2_pathway) = c2_pathway[,1]
  c2_list = lapply(c2_list, function(x) x[-(1:2)])
  names(c2_list) = c2_pathway[,1]
  c2_len = lapply(c2_list, length)
  gene2path = data.frame(Gene = unlist(c2_list), PathwayID = rep(names(c2_list), c2_len),
                         stringsAsFactors = FALSE)
  limit_pathways = names(c2_list)[c2_len>=limit[1] & c2_len<=limit[2]]
  gene2path = gene2path[gene2path$PathwayID%in%limit_pathways, ]
  gene2path$PathwayName = c2_pathway[gene2path$PathwayID, 2]
  return(gene2path)
}


#' Write GMT file
#'
#' write data frame to a gmt file
#'
#' @docType methods
#' @name ReadGMT
#' @rdname ReadGMT
#'
#' @param gene2path A data frame. The columns should be Gene, Pathway ID, and Pathway Name.
#' @param gmtfile Path to gmt file.
#'
#' @return Output gmt file to local folder.
#'
#' @author Wubing Zhang
#' @examples
#' gene2path = gsGetter(type = "Complex")
#' # writeGMT(gene2path, "Protein_complex.gmt")
#'
#' @export
#'
writeGMT <- function(gene2path, gmtfile){
  genelists = lapply(unique(gene2path[,2]), function(x)
    paste(gene2path[gene2path[,2]==x, 1], collapse = "\t"))
  gmt = cbind(unique(gene2path[,2]), gene2path[!duplicated(gene2path[,2]),3], unlist(genelists))
  write.table(gmt, gmtfile, sep = "\t", row.names = FALSE,
              col.names = FALSE, quote = FALSE)
}
