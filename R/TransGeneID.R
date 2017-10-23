TransGeneID <- function(genes, fromType="SYMBOL", toType="ENTREZID", organism="hsa"){
  if(is.null(genes) | length(genes)==0){
    return(c())
  }
  genes = as.character(genes)
  tmpGene = toupper(genes)
  #=====================
  GeneID_Convert = getOrg(organism, onlyLib = FALSE)$Symbol_Entrez
  idx = duplicated(GeneID_Convert[,fromType])
  convert = GeneID_Convert[!idx,toType]
  names(convert) = GeneID_Convert[!idx,fromType]
  gene_after=convert[tmpGene]
  names(gene_after)=genes

  return(gene_after)
}
