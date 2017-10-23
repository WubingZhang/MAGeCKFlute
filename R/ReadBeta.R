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
    stop("The parameter gene_summary is not permitted!")
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
  dd1$Control = rowMeans(dd[,ctrlName,drop= FALSE])
  dd1$Treatment = rowMeans(dd[,treatName,drop= FALSE])
  dd1$ENTREZID = TransGeneID(dd$Gene, fromType = "SYMBOL",
                             toType = "ENTREZID", organism = organism)[dd$Gene]
  dd1=as.data.frame(dd1, stringsAsFactors= FALSE)

  ##==============Remove NAs=============================================
  idx=is.na(dd1$ENTREZID)
  dd1=dd1[!idx,]

  return(dd1)
}


