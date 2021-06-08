#' Load processed Depmap data
#'
#' @docType methods
#' @name loadDepmap
#' @rdname loadDepmap
#'
#' @author Wubing Zhang
#'
#' @return A list including two elements, one is the Depmap CRISPR data, and the
#' other is the sample annotation data.
#'
#' @examples
#' \dontrun{
#'   depmapDat = LoadDepmap()
#' }
#'
#' @export

LoadDepmap <- function(){
  ## Load Depmap data
  depmap_rds = file.path(system.file("extdata", package = "MAGeCKFlute"), "Depmap.rds")
  if(file.exists(depmap_rds)){
    Depmap = readRDS(depmap_rds)
  }else{
    Depmap = t(read.csv("https://ndownloader.figshare.com/files/25494359", header = TRUE,
                        row.names = 1, stringsAsFactors = FALSE, check.names = FALSE))
    rownames(Depmap) = gsub(" .*", "", rownames(Depmap))
    saveRDS(Depmap, depmap_rds)
  }
  meta_rds = file.path(system.file("extdata", package = "MAGeCKFlute"), "Depmap_sample_info.rds")
  if(file.exists(meta_rds)){
    sampleinfo = readRDS(meta_rds)
  }else{
    sampleinfo = read.csv("https://ndownloader.figshare.com/files/25494443",
                          row.names = 1, header = TRUE, stringsAsFactors = FALSE)
    saveRDS(sampleinfo, meta_rds)
  }
  return(list(Depmap = Depmap, sampleinfo = sampleinfo))
}

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
#' file1 = file.path(system.file("extdata", package = "MAGeCKFlute"),
#' "testdata/rra.gene_summary.txt")
#' gdata = ReadRRA(file1)
#' head(gdata)
#' \dontrun{
#'   gdata = IncorporateDepmap(gdata)
#'   head(gdata)
#' }
#' @export
#'
IncorporateDepmap <- function(dd, symbol = "id",
                              cell_lines = NA, lineages = "All",
                              na.rm = FALSE){
  depmapDat = LoadDepmap()
  Depmap = depmapDat$Depmap
  sampleinfo = depmapDat$sampleinfo

  sampleinfo = sampleinfo[colnames(Depmap), ]
  idx1 = which(sampleinfo$lineage%in%lineages)
  idx2 = which(sampleinfo$stripped_cell_line_name%in%cell_lines)
  idx3 = which(sampleinfo$CCLE.Name%in%cell_lines)
  idx4 = which(sampleinfo$alias%in%cell_lines)
  idx2 = c(idx2, idx3, idx4)
  genes = as.character(dd[, symbol])
  Depmap = as.data.frame(Depmap, stringsAsFactors = FALSE)
  if(length(idx2)>0){
    dd = cbind(dd, Depmap = rowMeans(Depmap[genes, idx2, drop=FALSE], na.rm = TRUE))
  }else if(length(idx1)>0){
    dd = cbind(dd, Depmap = rowMeans(Depmap[genes, idx1, drop=FALSE], na.rm = TRUE))
  }else{
    dd = cbind(dd, Depmap = rowMeans(Depmap[genes, ], na.rm = TRUE))
  }
  if(na.rm) dd = na.omit(dd)
  return(dd)
}


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
#' file1 = file.path(system.file("extdata", package = "MAGeCKFlute"),
#' "testdata/rra.gene_summary.txt")
#' gdata = ReadRRA(file1)
#' \dontrun{
#'   rra.omit = OmitCommonEssential(gdata)
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
    Depmap_19Q3 = t(read.csv("https://ndownloader.figshare.com/files/24613292", header = TRUE,
                             row.names = 1, stringsAsFactors = FALSE, check.names = FALSE))
    rownames(Depmap_19Q3) = gsub(" .*", "", rownames(Depmap_19Q3))
    saveRDS(Depmap_19Q3, depmap_rds)
  }
  meta_rds = file.path(system.file("extdata", package = "MAGeCKFlute"), "Depmap_sample_info.rds")
  if(file.exists(meta_rds)){
    sampleinfo = readRDS(meta_rds)
  }else{
    sampleinfo = read.csv("https://ndownloader.figshare.com/files/24613394",
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
  depmapDat = LoadDepmap()
  Depmap = depmapDat$Depmap
  sampleinfo = depmapDat$sampleinfo
  if(!"all" %in% tolower(lineages)){
    idx = sampleinfo$lineage%in%tolower(lineages)
    idx = colnames(Depmap)%in%rownames(sampleinfo)[idx]
    if(sum(idx)>5){
      Depmap = Depmap[, idx]
    }else{ warning("Less than 5 cell lines are avaible, so ignore lineage setting.")}
  }
  idx = rowSums(Depmap<dependency, na.rm = TRUE)>0.6*ncol(Depmap)
  lethal_genes = rownames(Depmap)[idx]
  dd = dd[!(dd[,symbol] %in% lethal_genes), ]
  return(dd)
}
