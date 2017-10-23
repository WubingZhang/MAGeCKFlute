#===Distribution of beta scores======================================
CellCycleFit <- function(beta,ctrlname="Control",treatname="Treatment",
                         main=NULL,ylab="Beta score",filename=NULL){
  loginfo(paste("Cell cycle fitting for", main, ylab))
  dd2 = beta[,c("Gene", ctrlname, treatname)]
  dd2 = melt(dd2,id="Gene")
  dd2$x = rep(rowMeans(beta[,ctrlname,drop=FALSE]),ncol(dd2)-1)

  p=ggplot(dd2,aes(x,value,color=variable,group=variable))
  p=p+geom_point(alpha=4/10,size=0.8)
  p=p+geom_smooth(method='lm',se=FALSE)
  p=p+theme_bw(12)+theme(plot.title = element_text(hjust = 0.5,size=12))
  p=p+labs(x="Control",y=ylab,title=main,color=NULL)#+ggtitle("Normalization with")
  p=p+theme(legend.justification = c(0, 1), legend.position = c(0.01, 0.99))

  if(!is.null(filename)){
    ggsave(plot=p,filename=filename,units = "in",width=300/100,height =270/100 )
  }
  return(p)
}


