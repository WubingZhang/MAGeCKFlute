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
#' @import depmap ExperimentHub
#' @export

LoadDepmap <- function(){
  requireNamespace("depmap", quietly=TRUE) || stop("need depmap package")
  requireNamespace("ExperimentHub", quietly=TRUE) || stop("need ExperimentHub package")
  eh <- ExperimentHub()
  crispr <- eh[["EH2261"]]
  metadata <- eh[["EH2266"]]
  metadata = as.data.frame(metadata)
  rownames(metadata) = metadata$depmap_id
  metadata = metadata[,-1]
  crispr = matrix(crispr$dependency, nrow = length(unique(crispr$gene_name)),
                  ncol = length(unique(crispr$depmap_id)),
                  dimnames = list(rownames = unique(crispr$gene_name),
                                  colnames = unique(crispr$depmap_id)))
  crispr = as.data.frame(crispr)
  # ## Load Depmap data
  # depmap_rds = file.path(system.file("extdata", package = "MAGeCKFlute"), "Depmap.rds")
  # if(file.exists(depmap_rds)){
  #   Depmap = readRDS(depmap_rds)
  # }else{
  #   Depmap = t(read.csv("https://ndownloader.figshare.com/files/25494359", header = TRUE,
  #                       row.names = 1, stringsAsFactors = FALSE, check.names = FALSE))
  #   rownames(Depmap) = gsub(" .*", "", rownames(Depmap))
  #   saveRDS(Depmap, depmap_rds)
  # }
  # meta_rds = file.path(system.file("extdata", package = "MAGeCKFlute"), "Depmap_sample_info.rds")
  # if(file.exists(meta_rds)){
  #   sampleinfo = readRDS(meta_rds)
  # }else{
  #   sampleinfo = read.csv("https://ndownloader.figshare.com/files/25494443",
  #                         row.names = 1, header = TRUE, stringsAsFactors = FALSE)
  #   saveRDS(sampleinfo, meta_rds)
  # }
  return(list(Depmap = crispr, sampleinfo = metadata))
}

#' Incorporate Depmap screen into analysis
#'
#' @docType methods
#' @name IncorporateDepmap
#' @rdname IncorporateDepmap
#'
#' @param dd A data frame.
#' @param symbol A character, specifying the column name of gene symbols in the data frame.
#' @param cell_lines A character vector, specifying the cell lines for incorporation.
#' @param lineages A character vector, specifying the cancer types for incorporation.
#' @param na.rm Boolean, indicating whether removing NAs from the results.
#' @author Wubing Zhang
#'
#' @return A data frame with Depmap column (average CERES scores across selected cell lines) attached.
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
  sampleinfo$primary_disease = tolower(sampleinfo$primary_disease)
  sampleinfo$subtype_disease = tolower(sampleinfo$subtype_disease)
  sampleinfo$cell_line = tolower(sampleinfo$cell_line)
  sampleinfo = sampleinfo[colnames(Depmap), ]

  idx0 = which(sampleinfo$primary_disease %in% tolower(lineages))
  idx1 = which(sampleinfo$subtype_disease %in% tolower(lineages))
  if("all" %in% tolower(lineages)){
    idx1 = 1:ncol(Depmap)
  }else{
    idx1 = unique(c(idx0, idx1))
  }
  idx2 = which(gsub("_.*", "", sampleinfo$cell_line) %in% tolower(cell_lines))
  idx3 = which(sampleinfo$cell_line%in%tolower(cell_lines))
  idx4 = which(sampleinfo$aliases%in%cell_lines)
  idx2 = unique(c(idx2, idx3, idx4))
  genes = as.character(dd[, symbol])
  if(length(idx2)>0){
    dd = cbind(dd, Depmap = rowMeans(Depmap[genes, idx2, drop=FALSE], na.rm = TRUE))
  }else if(length(idx1)>0){
    dd = cbind(dd, Depmap = rowMeans(Depmap[genes, idx1, drop=FALSE], na.rm = TRUE))
  }else{
    warning("No cell line is selected. Using the average score across all cell types")
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
  depmapDat = LoadDepmap()
  Depmap = depmapDat$Depmap
  sampleinfo = depmapDat$sampleinfo
  sampleinfo$primary_disease = tolower(sampleinfo$primary_disease)
  sampleinfo$subtype_disease = tolower(sampleinfo$subtype_disease)
  sampleinfo$cell_line = tolower(sampleinfo$cell_line)
  sampleinfo = sampleinfo[colnames(Depmap), ]

  if(!"all" %in% tolower(lineages)){
    idx0 = which(sampleinfo$primary_disease %in% tolower(lineages))
    idx1 = which(sampleinfo$subtype_disease %in% tolower(lineages))
    idx = unique(c(idx0, idx1))
    if(length(idx)>5){
      Depmap = Depmap[, idx]
    }else{ warning("Less than 5 cell lines are avaible, so ignore lineage setting.")}
  }
  ## Explore the relationship
  genes = intersect(rownames(dd), rownames(Depmap))
  if(length(genes)<10) stop("Invalid gene symbols.")
  if(method %in% c("pearson", "spearman", "kendall")){
    similarity = apply(Depmap[genes,], 2, function(x){
      tmp = cor.test(x, dd[genes, score], method = method, na.action=na.omit)
      c(tmp$estimate, tmp$p.value)
    })
  }else{
    stop("Invalid distance measure!!!")
  }
  similarity = as.data.frame(t(similarity))
  colnames(similarity) = c("estimate", "p.value")
  rownames(similarity) = sampleinfo[colnames(Depmap), 1]
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
#' @param lineages A character vector, specifying the lineages for selecting essential genes.
#' @param cell_lines A character vector, specifying cell lines for selecting essential genes.
#' @param dependency A numeric, specifying the threshold for selecting essential genes.
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
                                cell_lines = NULL,
                                dependency = -0.5){
  depmapDat = LoadDepmap()
  Depmap = depmapDat$Depmap
  sampleinfo = depmapDat$sampleinfo
  sampleinfo$primary_disease = tolower(sampleinfo$primary_disease)
  sampleinfo$subtype_disease = tolower(sampleinfo$subtype_disease)
  sampleinfo$cell_line = tolower(sampleinfo$cell_line)
  sampleinfo = sampleinfo[colnames(Depmap), ]

  idx0 = which(sampleinfo$primary_disease %in% tolower(lineages))
  idx1 = which(sampleinfo$subtype_disease %in% tolower(lineages))
  if("all" %in% tolower(lineages)){
    idx1 = 1:ncol(Depmap)
  }else{
    idx1 = unique(c(idx0, idx1))
  }
  idx2 = which(gsub("_.*", "", sampleinfo$cell_line) %in% tolower(cell_lines))
  idx3 = which(sampleinfo$cell_line%in%tolower(cell_lines))
  idx4 = which(sampleinfo$aliases%in%cell_lines)
  idx2 = unique(c(idx2, idx3, idx4))
  if(length(idx2)>0){
    Depmap = Depmap[, idx2, drop=FALSE]
  }else if(length(idx1)>0){
    Depmap = Depmap[, idx1, drop=FALSE]
  }else{
    warning("No cell line is selected. Using all cell types instead")
  }

  idx = rowSums(Depmap<dependency, na.rm = TRUE)>0.6*ncol(Depmap)
  lethal_genes = rownames(Depmap)[idx]
  dd = dd[!(dd[,symbol] %in% lethal_genes), ]
  return(dd)
}
