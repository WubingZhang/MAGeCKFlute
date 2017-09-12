
#===Distribution of beta scores======================================
Density.plot <- function(beta, ctrlname="Control",treatname="Treatment",
                         main=NULL,xlab="Beta Score",filename=NULL){
  dd1 = beta
  loginfo(paste("Density plot for", main, xlab, "..."))
  dd1 = dd1[,c("Gene", ctrlname, treatname)]
  dd1 = melt(dd1[,1:3],id="Gene")
  dd1$variable = gsub(".beta","",dd1$variable)
  p=ggplot(data=dd1,aes(x=value,color=variable,group=variable))
  p=p+geom_density()
  # p=p+facet_wrap(~variable,nrow=1)
  p=p+scale_color_npg()+theme_bw(14)+theme(plot.title = element_text(hjust = 0.5,size=12))
  p=p+labs(color=NULL)
  p=p+theme(legend.justification = c(0, 1), legend.position = c(0.01, 0.99))
  # p=p+theme(legend.text = element_text(size=8))
  p=p+labs(x=xlab, y="Density", title=main)

  if(!is.null(filename)){
    ggsave(plot=p,filename=filename,units = "in",width=300/100,height =270/100 )
  }
  return(p)
}
