#' Read gene beta scores
#'
#' Read gene beta scores from file or data frame
#'
#' @docType methods
#' @name ReadBeta
#' @rdname ReadBeta
#' @aliases readbeta
#'
#' @param gene_summary a file path or a data frame, data frame, which has columns of 'Gene',
#' \code{ctrlname}.beta and \code{treatname}.beta.
#' @param ctrlName character vector, specifying the name of control sample.
#' @param treatName character vector, specifying the name of treatment sample.
#' @param organism character, KEGG species code, or the common species name, used to determine
#' the gene annotation package. For all potential values check: data(bods); bods. Default org="hsa",
#' and can also be "human" (case insensitive).
#'
#' @return a data frame including four columns, named "Gene", "Control", "Treatment", and "ENTREZID".
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
#' ctrlName = c("D7_R1", "D7_R2")
#' treatName = c("PLX7_R1", "PLX7_R2")
#' #Read beta score from gene summary table
#' dd = ReadBeta(MLE_Data, ctrlName = ctrlName, treatName = treatName, organism="hsa")
#' head(dd)
#'
#' @import pathview
#'
#' @export


#===read gene summary file=============================================
ReadBeta <- function(gene_summary, ctrlName="Control",
                     treatName="Treatment", organism='hsa'){
  loginfo("Read gene summary file ...")

  #=========If gene_summary is a path or a data frame====================
  if(class(gene_summary)=="character" && file.exists(gene_summary)){
    dd=read.table(file=gene_summary,header= TRUE, stringsAsFactors = FALSE)
  }else if(class(gene_summary)=="data.frame" &&
           all(c("Gene", paste(c(ctrlName, treatName),"beta",sep ="."))
               %in%colnames(gene_summary))){
    dd = gene_summary
  }else{
    stop("gene_summary is invalid!")
  }

  #=========Remove non-target control sgRNA==============================
  idx=grepl("^CTRL",dd$Gene,ignore.case = TRUE)
  dd=dd[!idx,]
  idx=grepl("beta",names(dd))
  idx[1]= TRUE
  dd=dd[,idx]
  names(dd)=gsub(".beta","",names(dd))

  #==============Deal with replicates and convert geneid==================
  dd1 = list()
  dd1$Gene = dd$Gene
  dd1$Control = rowMeans(dd[,ctrlName, drop = FALSE])
  dd1$Treatment = rowMeans(dd[,treatName, drop = FALSE])
  dd1$ENTREZID = suppressMessages(id2eg(dd$Gene, "SYMBOL", org = organism)[, "ENTREZID"])
  dd1=as.data.frame(dd1, stringsAsFactors= FALSE)

  ##==============Remove NAs=============================================
  idx=is.na(dd1$ENTREZID)
  dd1=dd1[!idx,]
  return(dd1)
}


