#' Gene ID conversion between ENTREZID and SYMBOL
#'
#' Gene ID conversion between ENTREZID and SYMBOL
#'
#' @docType methods
#' @name TransGeneID
#' @rdname TransGeneID
#' @aliases transGeneID
#'
#' @param genes A character vector, input genes to be converted.
#' @param fromType The input ID type, one of "Symbol" (default), "Entrez",
#' "Ensembl", "FULLNAME" and "TYPE". you can e.g. do: mart = useMart("ensembl", dataset="hsapiens_gene_ensembl"),
#' followed by listAttributes(mart) to get valid attribute names.
#' @param toType The output ID type, one of "Symbol", "Entrez" (default),
#' "Ensembl", "FULLNAME" and "TYPE".
#' @param organism A character, specifying organism, one of "hsa"(or 'Human'),
#' "mmu"(or 'Mouse'), "bta", "cfa", "ptr", "rno", and "ssc".
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
#' data(MLE_Data)
#' TransGeneID(MLE_Data$Gene[1:10], organism="hsa")
#' TransGeneID(MLE_Data$Gene[1:10], organism="hsa", useBiomart = TRUE,
#'             fromType = "hgnc_symbol", toType = "entrezgene")
#' @import biomaRt
#' @export


TransGeneID <- function(genes, fromType="Symbol", toType="Entrez",
                        organism="hsa", useBiomart = FALSE,
                        ensemblHost = "www.ensembl.org"){
  requireNamespace("biomaRt")
  if(is.null(genes) | length(genes)==0){
    return(c())
  }
  genes = as.character(genes)
  tmpGene = genes
  #=============Read annotation file========
  orgInfo <- getOrg(organism)
  GeneID_Convert = orgInfo$Symbol_Entrez
  if(useBiomart){
    # listEnsemblArchives()
    # listMarts()
    mart = useMart('ensembl')
    datasets = listDatasets(mart)
    data(bods)
    tmp = bods[bods[,3] == orgInfo$org, 2]
    ds = datasets[grepl(tmp, datasets$description, ignore.case = TRUE), 1][1]
    ensembl <- useMart(host = ensemblHost, biomart = 'ENSEMBL_MART_ENSEMBL', dataset = ds)
    GeneID_Convert = getBM(attributes=c(fromType, toType), mart = ensembl)
  }
  ##======Convert==========
  fromType = tolower(fromType)
  toType = tolower(toType)
  idx = duplicated(GeneID_Convert[, fromType])
  convert = GeneID_Convert[!idx, toType]
  names(convert) = GeneID_Convert[!idx,fromType]
  gene_after=convert[tmpGene]
  names(gene_after)=genes
  return(gene_after)
}
