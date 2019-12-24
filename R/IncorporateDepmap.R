#' Incorporate Depmap screen into analysis
#'
#' @docType methods
#' @name IncorporateDepmap
#' @rdname IncorporateDepmap
#'
#' @param dd A data frame.
#' @param symbol A character, specifying the column name of gene symbols in the data frame.
#' @param cell_lines A character vector, specifying the cell lines in Depmap to be considered.
#' @param lineages A character vector, specifying the lineages in Depmap to be considered.
#' @param na.rm Boolean, indicating whether removing NAs from the results.
#' @author Wubing Zhang
#'
#' @return A data frame with Depmap column attached.
#'
#' @examples
#' dd.rra = ReadRRA(rra.gene_summary)
#' depmap_similarity = ResembleDepmap(dd.rra)
#' dd.rra = IncorporateDepmap(dd.rra, cell_lines=rownames(depmap_similarity)[1:3])
#' head(dd.rra)
#' @export
#'
IncorporateDepmap <- function(dd, symbol = "id", cell_lines = NA, lineages = "All", na.rm = TRUE){
  depmap_rds = file.path(system.file("extdata", package = "MAGeCKFlute"), "Depmap_19Q3.rds")
  if(file.exists(depmap_rds)){
    Depmap_19Q3 = readRDS(depmap_rds)
  }else{
    locfname <- file.path(system.file("extdata", package = "MAGeCKFlute"), "Achilles_gene_effect.csv")
    download.file("https://ndownloader.figshare.com/files/20234073", locfname, quiet = FALSE)
    Depmap_19Q3 = t(read.csv(locfname, header = TRUE, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE))
    rownames(Depmap_19Q3) = gsub(" .*", "", rownames(Depmap_19Q3))
    saveRDS(Depmap_19Q3, depmap_rds)
  }
  sampleinfo = readRDS(file.path(system.file("extdata", package = "MAGeCKFlute"), "Depmap_sample_info.rds"))
  sampleinfo = sampleinfo[colnames(Depmap_19Q3), ]
  idx1 = sampleinfo$lineage%in%lineages
  idx2 = sampleinfo$stripped_cell_line_name%in%cell_lines |
    sampleinfo$CCLE.Name%in%cell_lines |
    sampleinfo$alias%in%cell_lines
  genes = as.character(dd[, symbol])
  Depmap_19Q3 = as.data.frame(Depmap_19Q3, stringsAsFactors = FALSE)
  if(sum(idx2)>0){
    dd = cbind(dd, Depmap = rowMeans(Depmap_19Q3[genes, idx2], na.rm = TRUE))
  }else if(sum(idx1)>0){
    dd = cbind(dd, Depmap = rowMeans(Depmap_19Q3[genes, idx1], na.rm = TRUE))
  }else{
    dd = cbind(dd, Depmap = rowMeans(Depmap_19Q3[genes, ], na.rm = TRUE))
  }
  if(na.rm) dd = na.omit(dd)
  return(dd)
}
