
#===Distribution of beta scores======================================
Violin.plot <- function(beta, ctrlname="Control",treatname="Treatment",
                        main=NULL,ylab="Beta Score",filename=NULL){

  dd1=beta
  loginfo(paste("Violin plot for", main, ylab, "..."))
  dd1 = dd1[,c("Gene", ctrlname, treatname)]
  dd1 = melt(dd1[,1:3],id="Gene")
  dd1$variable = gsub(".beta","",dd1$variable)
  p=ggplot(data=dd1,aes(x=variable,y=value,color=variable))
  p=p+geom_violin()+geom_boxplot(width=.1, outlier.colour=NA)
  #p=p+ylim(-1.5,1)
  p=p+scale_color_npg()+theme_bw(14)+theme(plot.title = element_text(hjust = 0.5,size=12))
  p=p+theme(legend.position = "none")
  p=p+theme(axis.text.x = element_text(angle = 30, hjust = 1, size=8))
  p=p+labs(x="Sample",y=ylab,title=main)

  if(!is.null(filename)){
    ggsave(plot=p,filename=filename,units = "in",width=300/100,height =270/100 )
  }
  return(p)
}
