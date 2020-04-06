#' Compute the similarity between customized CRISPR screen with Depmap screens
#'
#' @docType methods
#' @name ResembleDepmap
#' @rdname ResembleDepmap
#'
#' @param dd A data frame.
#' @param symbol A character, specifying the column name of gene symbols in the data frame.
#' @param score A character, specifying the column name of gene essentiality score in the data frame.
#' @param lineages A character vector, specifying the lineages used for common essential gene selection.
#' @param method A character, indicating which correlation coefficient is to be used for the test. One of "pearson", "kendall", or "spearman".
#'
#' @author Wubing Zhang
#'
#' @return A data frame with correlation and test p.value.
#'
#' @examples
#' dd.rra = ReadRRA(rra.gene_summary)
#' \dontrun{
#'   rra.omit = OmitCommonEssential(dd.rra)
#'   depmap_similarity = ResembleDepmap(rra.omit)
#'   head(depmap_similarity)
#' }
#' @export

ResembleDepmap <- function(dd, symbol = "id", score = "Score", lineages = "All",
                           method = c("pearson", "spearman", "kendall")[1]){
  dd = dd[!duplicated(dd[, symbol]), ]
  rownames(dd) = dd[, symbol]

  ## Load Depmap data
  depmap_rds = file.path(system.file("extdata", package = "MAGeCKFlute"), "Depmap_19Q3.rds")
  if(file.exists(depmap_rds)){
    Depmap_19Q3 = readRDS(depmap_rds)
  }else{
    Depmap_19Q3 = t(read.csv("https://ndownloader.figshare.com/files/20234073", header = TRUE,
                             row.names = 1, stringsAsFactors = FALSE, check.names = FALSE))
    rownames(Depmap_19Q3) = gsub(" .*", "", rownames(Depmap_19Q3))
    saveRDS(Depmap_19Q3, depmap_rds)
  }
  meta_rds = file.path(system.file("extdata", package = "MAGeCKFlute"), "Depmap_sample_info.rds")
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
  ## Explore the relationship
  genes = intersect(rownames(dd), rownames(Depmap_19Q3))
  if(length(genes)<10) stop("Invalid gene symbols.")
  if(method %in% c("pearson", "spearman", "kendall")){
    similarity = apply(Depmap_19Q3[genes,], 2, function(x){
      tmp = cor.test(x, dd[genes, score], method = method, na.action=na.omit)
      c(tmp$estimate, tmp$p.value)
    })
  }else{
    stop("Invalid distance measure!!!")
  }
  similarity = as.data.frame(t(similarity))
  colnames(similarity) = c("estimate", "p.value")
  rownames(similarity) = sampleinfo[colnames(Depmap_19Q3), 1]
  similarity = similarity[order(-similarity$estimate), ]
  return(similarity)
}
