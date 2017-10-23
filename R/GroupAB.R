GroupAB <- function(beta, ctrlname="Control", treatname="Treatment",
                    scale_cutoff=1, filename=NULL){
  loginfo("Select AB group genes ...")
  dd1 = beta[,c("Gene",'ENTREZID',treatname, ctrlname)]
  dd1$diff=rowMeans(dd1[,treatname,drop= FALSE])-rowMeans(dd1[,ctrlname,drop= FALSE])
  cutoff=Cutoff_Calling(dd1$diff,scale=scale_cutoff)
  dd1$group="no"
  dd1$group[dd1$diff>cutoff]="up"
  dd1$group[dd1$diff<(-cutoff)]="down"

  if(!is.null(filename)){
    write.table(dd1, filename, sep="\t", row.names = FALSE,col.names = TRUE,quote= FALSE)
  }
  return(dd1)
}
