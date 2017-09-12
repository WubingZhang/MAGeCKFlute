#===read gene summary file=====================================
ReadRRA <- function(gene_summary){
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

  # tmp <- as.data.frame(org.Hs.egALIAS2EG)
  data("org.Hs.egALIAS2EG")
  idx = which(tmp$alias_symbol %in% dd$id)
  tmp = tmp[idx,]
  tmp1=suppressMessages(as.data.frame(bitr(tmp$gene_id, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")))
  tmp=merge(tmp,tmp1,by.x="gene_id",by.y="ENTREZID")
  tmp=tmp[,2:3]
  idx=duplicated(as.character(tmp$alias_symbol))
  tmp=tmp[!idx,]
  idx=duplicated(as.character(tmp$SYMBOL))
  tmp=tmp[!idx,]
  rownames(tmp)=tmp$alias_symbol

  dd$Official=tmp[as.character(dd$id),"SYMBOL"]
  idx=is.na(dd$Official)
  dd=dd[!idx,]
  tmp=suppressMessages(as.data.frame(bitr(dd$Official, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")))
  dd=merge(dd,tmp,by.x="Official",by.y="SYMBOL")
  idx=duplicated(dd$Official)
  dd=dd[!idx,]

  dd=dd[,c(1,3:5)]
  colnames(dd) = c("Official", "neg.fdr", "pos.fdr", "ENTREZID")
  return(dd)
}
