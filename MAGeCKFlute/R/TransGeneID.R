TransGeneID <- function(genes, fromType="SYMBOL", toType="ENTREZID", organism="hsa"){
  if(is.null(genes) | length(genes)==0){
    return(NULL)
  }
  genes = as.character(genes)
  #=====================
  if(organism=="hsa"){
    idx = duplicated(GeneID_Convert_hsa[,fromType])
    convert = GeneID_Convert_hsa[!idx,toType]
    names(convert) = GeneID_Convert_hsa[!idx,fromType]
    gene_after=convert[genes]
  }
  if(organism=="mmu"){
    tmpGenes = toupper(genes)
    idx = duplicated(GeneID_Convert_mmu[,fromType])
    convert = GeneID_Convert_mmu[!idx,toType]
    names(convert) = GeneID_Convert_mmu[!idx,fromType]
    gene_after=convert[tmpGenes]
    names(gene_after)=genes
  }
  return(gene_after)
}
