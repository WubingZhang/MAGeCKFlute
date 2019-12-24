#' Omit common essential genes based on depmap data
#'
#' @docType methods
#' @name OmitCommonEssential
#' @rdname OmitCommonEssential
#'
#' @param dd A data frame.
#' @param symbol A character, specifying the column name of gene symbols in the data frame.
#' @param lineages A character vector, specifying the lineages used for common essential gene selection.
#' @param dependency A numeric, specifying the threshold for common essential gene selection.
#'
#' @author Wubing Zhang
#'
#' @return A data frame.
#'
#' @examples
#' dd.rra = ReadRRA(rra.gene_summary)
#' dim(dd.rra)
#' rra.omit = OmitCommonEssential(dd.rra)
#' dim(rra.omit)
#'
#' @export

OmitCommonEssential <- function(dd, symbol = "id", lineages = "All", dependency = -0.5){
  Depmap_19Q3 = readRDS(file.path(system.file("extdata", package = "MAGeCKFlute"), "Depmap_19Q3.rds"))
  if(!"all" %in% tolower(lineages)){
    sampleinfo = readRDS(file.path(system.file("extdata", package = "MAGeCKFlute"), "Depmap_sample_info.rds"))
    idx = sampleinfo$lineage%in%tolower(lineages)
    idx = colnames(Depmap_19Q3)%in%rownames(sampleinfo)[idx]
    if(sum(idx)>5){
      Depmap_19Q3 = Depmap_19Q3[, idx]
    }else{ warning("Less than 5 cell lines are avaible, so ignore lineage setting.")}
  }
  idx = rowSums(Depmap_19Q3<dependency, na.rm = TRUE)>0.7*ncol(Depmap_19Q3)
  lethal_genes = rownames(Depmap_19Q3)[idx]
  dd = dd[!(dd[,symbol] %in% lethal_genes), ]
  return(dd)
}
