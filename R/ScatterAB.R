
ScatterAB<-function(beta, ctrlname="Control",treatname="Treatment", scale_cutoff=1,
                     main=NULL,filename=NULL){
  beta$Control=rowMeans(beta[,ctrlname,drop=F])
  beta$Treatment=rowMeans(beta[,treatname,drop=F])
  intercept=Cutoff_Calling(beta$Treatment-beta$Control, scale=scale_cutoff)
  data=beta
  loginfo(paste("Scatter plot of", main, "Treat-Ctrl beta scores ..."))
  mycolour=c("no"="aliceblue",  "up"="#e41a1c","down"="#377eb8")
  xmin=min(data$Control)
  xmax=max(data$Control)
  ymin=min(data$Treatment)
  ymax=max(data$Treatment)
  #=========
  p=ggplot(data,aes(x=Control,y=Treatment,colour=group,fill=group))
  p=p+geom_point(position = "identity",shape=".",alpha=1/100,size = 0.01,show.legend = F)
  p=p+scale_color_manual(values=mycolour)
  p=p+geom_jitter(position = "jitter",show.legend = F)
  p=p+theme_bw(14)+theme(plot.title = element_text(hjust = 0.5,size=12))
  p=p+geom_abline(intercept = -intercept)
  p=p+geom_abline(intercept = +intercept)
  p=p+labs(x="Control beta.score",y="Treatment beta.score",title=main)
  p=p+annotate("text",color="#e41a1c",x=xmin, y=ymax,hjust = 0,
               label=paste("GroupA: ",as.character(dim(data[data$group=="up",])[1]),sep=""))
  p=p+annotate("text",color="#377eb8",x=xmax, y=ymin,hjust = 1,
               label=paste("GroupB: ",as.character(dim(data[data$group=="down",])[1]),sep=""))
  #============
  if(!is.null(filename)){
    ggsave(plot=p,filename=filename,units = "in",width=5,height =4 )
  }
  return(p)
}

