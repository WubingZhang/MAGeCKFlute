
ScatterAB<-function(beta, ctrlname="Control",treatname="Treatment",
                    intercept=Cutoff_Calling(beta[, treatname]-beta[, ctrlname]), slope=1,main=NULL,filename=NULL){
  data=beta
  loginfo(paste("Scatter plot of", main, "Treat-Ctrl beta scores ..."))
  mycolour=c("no"="aliceblue",  "up"="#e41a1c","down"="#377eb8")
  xmin=min(data[, ctrlname])
  xmax=max(data[, ctrlname])
  ymin=min(data[, treatname])
  ymax=max(data[, treatname])

  p=ggplot(data,aes(x=eval(parse(text = ctrlname)),y=eval(parse(text = treatname)),colour=group,fill=group))
  p=p+geom_point(position = "identity",shape=".",alpha=1/100,size = 0.01,show.legend = F)+scale_color_manual(values=mycolour)
  p=p+geom_jitter(position = "jitter",show.legend = F)
  p=p+theme_bw(14)+theme(plot.title = element_text(hjust = 0.5,size=12))
  p=p+geom_abline(slope=slope,intercept = -intercept)
  p=p+geom_abline(slope=slope,intercept = +intercept)
  p=p+labs(x="Control beta.score",y="Treatment beta.score",title=main)
  p=p+annotate("text",color="#e41a1c",label=paste("GroupA: ",as.character(dim(data[data$group=="up",])[1]),sep=""), x=xmin, y=ymax,hjust = 0)
  p=p+annotate("text",color="#377eb8",label=paste("GroupB: ",as.character(dim(data[data$group=="down",])[1]),sep=""), x=xmax, y=ymin,hjust = 1)
  if(!is.null(filename)){
    ggsave(plot=p,filename=filename,units = "in",width=5,height =4 )
  }
  return(p)
}

