
RankAB <- function(beta, genelist=c(), top=10, bottom=10, scale_cutoff=1, main=NULL,filename=NULL){
  loginfo(paste("Rank of", main, "Treat-Ctrl beta scores ..."))
  cutoff=Cutoff_Calling(beta$diff, scale=scale_cutoff)
  mycolour=c("no"="darkgray",  "up"="#e41a1c","down"="#377eb8")
  data = beta
  data$Rank=rank(data$diff)

  idx=(data$Rank<=bottom) | (data$Rank>(max(data$Rank)-top)) | (data$Gene %in% genelist)
  p=ggplot()
  p=p+geom_point(aes(x=diff,y=Rank,color=factor(data$group)),data=data,size=0.5)
  p=p+theme_bw(14)+theme(plot.title = element_text(hjust = 0.5,size=12))
  p=p+scale_color_manual(values=mycolour)
  p=p+geom_jitter(width = 0.1, height = 0.1)
  p=p+theme(panel.background = element_rect(fill='white', colour='black'))
  p=p+geom_vline(xintercept = 0,linetype = "dotted")+geom_vline(xintercept = c(cutoff,-cutoff),linetype = "dotted")
  p=p+geom_label_repel(aes(x=diff, y=Rank,fill=group,label = Gene),data=data[idx,],
                       fontface = 'bold', color = 'white',size = 2.5,force=3,
                       box.padding = unit(0.4, "lines"),segment.color = 'grey50',
                       point.padding = unit(0.3, "lines"))
  p=p+scale_fill_manual(values=mycolour)
  p=p+labs(x="Treatment-Control Beta Score",y="Rank",title=main)
  p=p+theme(legend.position="none")#+ylim(-1000,7000)
  if(!is.null(filename)){
    ggsave(plot=p,filename=filename,units = "in",width=5,height=4 )
  }
  return(p)
}
