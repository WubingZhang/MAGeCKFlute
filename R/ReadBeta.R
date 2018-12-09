#' Read gene beta scores
#'
#' Read gene beta scores from file or data frame
#'
#' @docType methods
#' @name ReadBeta
#' @rdname ReadBeta
#' @aliases readbeta
#'
#' @param gene_summary A file path or a data frame, data frame, which has columns of 'Gene' and '*|beta'.
#' @param keytype Type of gene id in `gene_summary`, which should be one of "Entrez" or "Symbol".
#' @param organism Character, KEGG species code, or the common species name, used to determine
#' the gene annotation package. For all potential values check: data(bods); bods. Default org="hsa",
#' and can also be "human" (case insensitive).
#'
#' @return A data frame, in which the first column is ENTREZID, and the later columns are beta score for each samples.
#'
#' @author Wubing Zhang
#'
#' @examples
#' data(mle.gene_summary)
#' dd = ReadBeta(mle.gene_summary, organism="hsa")
#' head(dd)
#'
#' @export
ReadBeta <- function(gene_summary, keytype = "Symbol", organism = 'hsa'){
  message(Sys.time(), " # Read gene summary file ...")

  #=========If gene_summary is a path or a data frame=====
  if(is.character(gene_summary) && file.exists(gene_summary)){
    dd=read.table(file=gene_summary,header= TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  }else if(is.data.frame(gene_summary) &&
          ("Gene"%in%names(gene_summary)) &&
          any(grepl(".beta", names(gene_summary)))){
    dd = gene_summary
  }else{
    stop("gene_summary is invalid!")
  }

  #=========Remove non-target control sgRNA==============================
  idx=grepl("^CTRL",dd$Gene,ignore.case = TRUE)
  dd=dd[!idx,]
  idx=grepl(".beta",names(dd))
  idx[1]= TRUE
  dd=dd[,idx]
  names(dd)=gsub(".beta","",names(dd))

  ##==============Remove NAs=============================================
  idx = is.na(dd$Gene) | duplicated(dd$Gene)
  dd = dd[!idx,]
  rownames(dd) = dd$Gene

  ##=============Gene ID conversion===============
  if(keytype == "Symbol")
    dd$EntrezID = TransGeneID(dd$Gene, "Symbol", "Entrez", organism = organism)
  else{
    dd$EntrezID = dd$Gene
    dd$Gene = TransGeneID(dd$Gene, "Entrez", "Symbol", organism = organism)
  }
  idx = is.na(dd$EntrezID) | duplicated(dd$EntrezID)
  dd = dd[!idx,]
  rownames(dd) = dd$EntrezID
  dd = dd[, c("Gene", "EntrezID", colnames(dd)[c(-1, -ncol(dd))])]
  return(dd)
}
