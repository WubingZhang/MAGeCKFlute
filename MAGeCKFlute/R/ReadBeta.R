#===read gene summary file=====================================
ReadBeta <- function(gene_summary, ctrlName="Control", treatName="Treatment"){
  loginfo("Read gene summary file ...")

  if(class(gene_summary)=="character" && file.exists(gene_summary)){
    dd=read.table(file=gene_summary,header=T)
  }else if(class(gene_summary)=="data.frame" &&
           all(c("Gene", paste(c(ctrlName, treatName),"beta",sep ="."))%in%colnames(gene_summary))){
    dd = gene_summary
  }else{
    stop("The parameter gene_summary is below standard!")
  }

  idx=grepl("Zhang_",dd$Gene,ignore.case = T)
  dd=dd[!idx,]
  idx=grepl("^CTRL",dd$Gene,ignore.case = T)
  dd=dd[!idx,]
  idx=grepl("beta",names(dd))
  idx[1]=T
  dd=dd[,idx]
  names(dd)=gsub(".beta","",names(dd))

  # tmp <- as.data.frame(org.Hs.egALIAS2EG)
  data("org.Hs.egALIAS2EG")
  idx = which(tmp$alias_symbol%in%dd$Gene)
  tmp = tmp[idx,]
  tmp1=suppressMessages(as.data.frame(bitr(tmp$gene_id, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")))
  tmp=merge(tmp,tmp1,by.x="gene_id",by.y="ENTREZID")
  tmp=tmp[,2:3]
  idx=duplicated(as.character(tmp$alias_symbol))
  tmp=tmp[!idx,]
  idx=duplicated(as.character(tmp$SYMBOL))
  tmp=tmp[!idx,]
  rownames(tmp)=tmp$alias_symbol

  dd$Offical=tmp[as.character(dd$Gene),"SYMBOL"]
  idx=is.na(dd$Offical)
  dd=dd[!idx,]
  tmp=suppressMessages(as.data.frame(bitr(dd$Offical, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")))
  dd=merge(dd,tmp,by.x="Offical",by.y="SYMBOL")
  idx=duplicated(dd$Offical)
  dd=dd[!idx,]
  dd=dd[,c("Offical",ctrlName,treatName,"ENTREZID")]

  dd1 = list()
  dd1$Gene = dd$Offical
  dd1$Control = rowMeans(dd[,ctrlName,drop=F])
  dd1$Treatment = rowMeans(dd[,treatName,drop=F])
  dd1$ENTREZID = dd$ENTREZID
  dd1=as.data.frame(dd1)

  return(dd1)
}


