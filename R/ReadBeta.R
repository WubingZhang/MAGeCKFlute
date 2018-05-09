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
#' @param organism Character, KEGG species code, or the common species name, used to determine
#' the gene annotation package. For all potential values check: data(bods); bods. Default org="hsa",
#' and can also be "human" (case insensitive).
#'
#' @return A data frame, in which the first column is ENTREZID, and the later columns are beta score for each samples.
#'
#' @author Wubing Zhang
#'
#' @note
#' The source can be found by typing \code{MAGeCKFlute:::ReadBeta}
#' or \code{getMethod("ReadBeta")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/ReadBeta.R}
#' Users should find it easy to customize this function.
#'
#' @examples
#' data(MLE_Data)
#' dd = ReadBeta(MLE_Data, organism="hsa")
#' head(dd)
#'
#' @import pathview
#'
#' @export


#===read gene summary file=============================================
ReadBeta <- function(gene_summary, organism='hsa'){
  loginfo("Read gene summary file ...")

  #=========If gene_summary is a path or a data frame====================
  if(class(gene_summary)=="character" && file.exists(gene_summary)){
    dd=read.table(file=gene_summary,header= TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  }else if(class(gene_summary)=="data.frame" &&
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
  dd$ENTREZID = TransGeneID(dd$Gene, "SYMBOL", "ENTREZID", organism = organism)
  #==============Deal with replicates and convert geneid==================
  # dd1 = list()
  # dd1$Gene = dd$Gene
  # dd1$Control = rowMeans(dd[,ctrlName, drop = FALSE])
  # dd1$Treatment = rowMeans(dd[,treatName, drop = FALSE])
  # dd1$ENTREZID = TransGeneID(dd$Gene, "SYMBOL", "ENTREZID", organism = organism)
  # dd1=as.data.frame(dd1, stringsAsFactors= FALSE)

  ##==============Remove NAs=============================================
  idx = is.na(dd$ENTREZID) | duplicated(dd$ENTREZID)
  dd = dd[!idx,]
  dd = dd[, c(1,ncol(dd), 2:(ncol(dd)-1))]
  rownames(dd) = dd$Gene
  dd = dd[, -1]
  return(dd)
}
