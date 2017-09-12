GroupAB <- function(beta, ctrlname="Control",treatname="Treatment", cutoff_scale=F,filename=NULL){
  loginfo("Select AB group genes ...")
  dd1 = beta
  dd1=dd1[,c("Gene",treatname, ctrlname, "ENTREZID")]
  dd1$diff=dd1[, treatname]-dd1[, ctrlname]
  cutoff=Cutoff_Calling(dd1$diff,scale=cutoff_scale)
  dd1$group="no"
  dd1$group[dd1$diff>cutoff]="up"
  dd1$group[dd1$diff<(-cutoff)]="down"

  if(!is.null(filename)){
    write.table(dd1, filename, sep="\t", row.names = F,col.names = T,quote=F)
  }
  return(list(result=dd1, cutoff=cutoff))
}
