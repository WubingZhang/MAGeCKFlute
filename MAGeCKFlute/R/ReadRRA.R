#===read gene summary file=====================================
ReadRRA <- function(gene_summary, organism="hsa"){
  loginfo("Read gene summary file ...")

  if(class(gene_summary)=="character" && file.exists(gene_summary)){
    dd=read.table(file=gene_summary,header=T)
  }else if(class(gene_summary)=="data.frame" &&
           all(c("id", "neg.fdr","pos.fdr")%in%colnames(gene_summary))){
    dd=gene_summary
  }else{
    stop("The parameter gene_summary is below standard!")
  }
  idx=grepl("Zhang_", dd$id, ignore.case = T)
  dd=dd[!idx,]
  idx=grepl("^CTRL", dd$id, ignore.case = T)
  dd=dd[!idx,]

  dd=dd[,c("id","neg.fdr","pos.fdr")]
  dd$ENTREZID = TransGeneID(dd$id, "SYMBOL", "ENTREZID", organism = organism)
  idx=is.na(dd$ENTREZID)
  dd=dd[!idx,]
  colnames(dd) = c("Official", "neg.fdr", "pos.fdr", "ENTREZID")

  return(dd)
}
