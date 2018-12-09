#' Gene ID conversion between ENTREZID and SYMBOL
#'
#' @docType methods
#' @name TransGeneID
#' @rdname TransGeneID
#' @aliases transGeneID
#'
#' @param genes A character vector, input genes to be converted.
#' @param fromType The input ID type, one of "Symbol" (default), "Entrez" and "Ensembl";
#' you can also input other valid attribute names for biomart.
#' @param toType The output ID type, one of "Symbol", "Entrez" (default), "Ensembl";
#' you can also input other valid attribute names for biomart.
#' @param organism One of "hsa"(or 'Human'), "mmu"(or 'Mouse'), "bta", "cfa", "ptr", "rno", and "ssc".
#' @param useBiomart Boolean, indicating whether use Biomart to do the transformation.
#' @param ensemblHost String, specifying ensembl host, you can use `listEnsemblArchives()`
#' to show all available Ensembl archives hosts.
#'
#' @return A character vector, named by unique input gene ids.
#'
#' @author Wubing Zhang
#'
#' @seealso \code{\link[pathview]{eg2id}}
#'
#' @examples
#' data(mle.gene_summary)
#' TransGeneID(mle.gene_summary$Gene[1:10], organism="hsa")
#' TransGeneID(mle.gene_summary$Gene[1:10], organism="hsa", useBiomart = TRUE)
#'
#' @import biomaRt
#' @export


TransGeneID <- function(genes, fromType="Symbol", toType="Entrez",
                        organism="hsa", useBiomart = FALSE,
                        ensemblHost = "www.ensembl.org"){
  requireNamespace("biomaRt")
  if(is.null(genes) | length(genes)==0){
    return(c())
  }
  data(bods)
  ridx=grep(tolower(organism), tolower(bods[,2:3])) %% nrow(bods)
  org = bods[ridx,3]
  genes = as.character(genes)
  tmpGene = genes
  fromType = tolower(fromType)
  toType = tolower(toType)
  #=============Read annotation file========
  if(useBiomart){
    datasets = paste0(c("hsapiens", "mmusculus", "btaurus", "cfamiliaris",
                        "ptroglodytes", "rnorvegicus", "sscrofa"), "_gene_ensembl")
    type = c("ensembl_gene_id", "entrezgene", "hgnc_symbol")
    if(org=="mmu") type[3] = "mgi_symbol"
    names(type) = c("ensembl", "entrez", "symbol")
    fromType = ifelse(fromType%in% names(type), type[fromType], fromType)
    toType = ifelse(toType%in% names(type), type[toType], toType)
    # listEnsemblArchives()
    # listMarts()
    ds = datasets[grepl(org, datasets)]
    ensembl <- useMart(host = ensemblHost, biomart = 'ENSEMBL_MART_ENSEMBL', dataset = ds)
    GeneID_Convert = getBM(attributes=c(fromType, toType), mart = ensembl)
  }else{
    GeneID_Convert <- getOrg(org)$Symbol_Entrez
  }
  ##======Convert==========
  idx = duplicated(GeneID_Convert[, fromType])
  convert = GeneID_Convert[!idx, toType]
  names(convert) = GeneID_Convert[!idx,fromType]
  gene_after=as.character(convert[tmpGene])
  names(gene_after)=genes
  return(gene_after)
}
