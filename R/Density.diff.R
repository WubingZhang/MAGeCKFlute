
#===Distribution of beta scores======================================

Density.diff <- function(beta, ctrlname="Control",treatname="Treatment",
                         main=NULL,filename=NULL){

  loginfo(paste("Density plot for", main, "treat-control beta scores..."))
  d=beta
  d$Diff=d[,treatname]-d[,ctrlname]
  d$r <- rnorm(length(d$Diff), mean=0, sd=sd(d$Diff)-0.01)
  p=ggplot(d,aes(x=Diff))
  p=p+geom_histogram(aes(y = ..density..),fill="gray90",binwidth=0.02)
  p=p+geom_density(colour="black")
  p=p+geom_density(aes(x=r,y=..density..),linetype="dashed",colour="red")
  p=p+geom_vline(xintercept = 0,linetype="dashed")
  p=p+theme_bw(14)+theme(plot.title = element_text(hjust = 0.5,size=12))
  p=p+labs(x="Treat-Control Beta Score",y="Density",title=main)#+ggtitle("Normalization with")

  if(!is.null(filename)){
    ggsave(plot=p,filename=filename,units = "in",width=300/100,height =270/100 )
  }
  return(p)
}
