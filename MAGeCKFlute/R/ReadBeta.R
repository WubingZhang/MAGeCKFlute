#===read gene summary file=====================================
ReadBeta <- function(gene_summary, ctrlName="Control", treatName="Treatment", organism='hsa'){
  loginfo("Read gene summary file ...")

  #=========
  #Decide if gene_summary is a path or a data frame
  if(class(gene_summary)=="character" && file.exists(gene_summary)){
    dd=read.table(file=gene_summary,header=T, stringsAsFactors = F)
  }else if(class(gene_summary)=="data.frame" &&
           all(c("Gene", paste(c(ctrlName, treatName),"beta",sep ="."))%in%colnames(gene_summary))){
    dd = gene_summary
  }else{
    stop("The parameter gene_summary is not permitted!")
  }

  #=========
  #Remove non-target control sgRNA
  idx=grepl("Zhang_",dd$Gene,ignore.case = T)
  dd=dd[!idx,]
  idx=grepl("^CTRL",dd$Gene,ignore.case = T)
  dd=dd[!idx,]
  idx=grepl("beta",names(dd))
  idx[1]=T
  dd=dd[,idx]
  names(dd)=gsub(".beta","",names(dd))

  dd1 = list()
  dd1$Gene = dd$Gene
  dd1$Control = rowMeans(dd[,ctrlName,drop=F])
  dd1$Treatment = rowMeans(dd[,treatName,drop=F])
  dd1$ENTREZID = TransGeneID(dd$Gene, fromType = "SYMBOL", toType = "ENTREZID",
                             organism = organism)[dd$Gene]
  dd1=as.data.frame(dd1, stringsAsFactors=F)

  idx=is.na(dd1$ENTREZID)
  dd1=dd1[!idx,]
  idx=duplicated(dd1$ENTREZID)
  dd1=dd1[!idx,]

  return(dd1)
}


