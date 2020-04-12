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
#' file1 = file.path(system.file("extdata", package = "MAGeCKFlute"),
#' "testdata/rra.gene_summary.txt")
#' gdata = ReadRRA(file1)
#' dim(gdata)
#' \dontrun{
#'   rra.omit = OmitCommonEssential(gdata)
#'   dim(rra.omit)
#' }
#'
#' @export

OmitCommonEssential <- function(dd, symbol = "id",
                                lineages = "All",
                                dependency = -0.5){
  ## Load Depmap data
  depmap_rds = file.path(system.file("extdata", package = "MAGeCKFlute"), "Depmap_19Q3.rds")
  if(file.exists(depmap_rds)){
    Depmap_19Q3 = readRDS(depmap_rds)
  }else{
    Depmap_19Q3 = t(read.csv("https://ndownloader.figshare.com/files/20234073",
                             header = TRUE, row.names = 1, stringsAsFactors = FALSE,
                             check.names = FALSE))
    rownames(Depmap_19Q3) = gsub(" .*", "", rownames(Depmap_19Q3))
    saveRDS(Depmap_19Q3, depmap_rds)
  }
  meta_rds = file.path(system.file("extdata", package = "MAGeCKFlute"),
                       "Depmap_sample_info.rds")
  if(file.exists(meta_rds)){
    sampleinfo = readRDS(meta_rds)
  }else{
    sampleinfo = read.csv("https://ndownloader.figshare.com/files/20274744",
                          row.names = 1, header = TRUE, stringsAsFactors = FALSE)
    saveRDS(sampleinfo, meta_rds)
  }
  if(!"all" %in% tolower(lineages)){
    idx = sampleinfo$lineage%in%tolower(lineages)
    idx = colnames(Depmap_19Q3)%in%rownames(sampleinfo)[idx]
    if(sum(idx)>5){
      Depmap_19Q3 = Depmap_19Q3[, idx]
    }else{ warning("Less than 5 cell lines are avaible, so ignore lineage setting.")}
  }
  idx = rowSums(Depmap_19Q3<dependency, na.rm = TRUE)>0.6*ncol(Depmap_19Q3)
  lethal_genes = rownames(Depmap_19Q3)[idx]
  dd = dd[!(dd[,symbol] %in% lethal_genes), ]
  return(dd)
}
