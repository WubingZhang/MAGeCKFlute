#Plot square
Square.plot<-function(beta, ctrlname="Control",treatname="Treatment",
                      main=NULL, filename=NULL, out.dir=".", scale_cutoff=1){
  loginfo(paste("Square plot for",main, "beta scores ..."))
  beta$group="Others"
  beta$Control = rowMeans(beta[, ctrlname, drop=F])
  beta$Treatment = rowMeans(beta[, treatname, drop=F])
  Control_cutoff = Cutoff_Calling(beta$Control, scale=scale_cutoff)
  drug_cutoff = Cutoff_Calling(beta$Treatment, scale=scale_cutoff)
  idx1 = -Control_cutoff<beta$Control
  idx2 = beta$Control<Control_cutoff
  idx3 = beta$Treatment>drug_cutoff
  idx=idx1&idx2&idx3
  beta$group[idx]="Group2" #proliferation
  idx3=beta$Treatment<(-drug_cutoff)
  idx=idx1&idx2&idx3
  beta$group[idx]="Group4" #anti-proliferation
  idx1=-drug_cutoff<beta$Treatment
  idx2=beta$Treatment<drug_cutoff
  idx3=beta$Control>Control_cutoff
  idx=idx1&idx2&idx3
  beta$group[idx]="Group3" #proliferation
  idx3=beta$Control<(-Control_cutoff)
  idx=idx1&idx2&idx3
  beta$group[idx]="Group1" #anti-proliferation
  #=========================
  x=beta$Control
  y=beta$Treatment
  slope=drug_cutoff/Control_cutoff
  intercept=0.2*Cutoff_Calling(y-x,scale=scale_cutoff)
  y_max=x*slope+intercept
  y_min=x*slope-intercept
  idx=y>y_max | y<y_min
  if(length(idx)>0)
    beta[!idx,"group"]="Others"
  #===================
  x_min=round(quantile(beta$Control[beta$group!="Others"],0.01),2)-0.01
  x_max=round(quantile(beta$Control[beta$group!="Others"],0.99),2)+0.01
  y_min=round(quantile(beta$Treatment[beta$group!="Others"],0.01),2)-0.01
  y_max=round(quantile(beta$Treatment[beta$group!="Others"],0.99),2)+0.01

  idx1=(x_min<=beta$Control & beta$Control<=x_max)
  idx2=(y_min<=beta$Treatment  & beta$Treatment<=y_max)
  idx=idx1&idx2

  gg=beta[idx,]
  # gg=beta
  gg$group=factor(gg$group,levels = c("Group1","Group2","Group3","Group4","Others"))
  #===============
  gg=gg[,c("Gene","Treatment","Control","group","ENTREZID")]
  mycolour=c("Others"="aliceblue",  "Group2"="#ff7f00", "Group3"="#ffff33",
             "Group4"="#984ea3", "Group1"="#4daf4a" )

  p=ggplot(gg,aes(x=Control,y=Treatment,colour=group,fill=group))
  p=p+geom_point(shape=".",alpha=1/1,size = 1)+scale_color_manual(values=mycolour)
  p=p+geom_jitter(size = 1)
  p=p+theme_bw(14)+theme(plot.title = element_text(hjust = 0.5,size=12))
  p=p+geom_vline(xintercept = Control_cutoff,linetype = "dotted")
  p=p+geom_vline(xintercept = -Control_cutoff,linetype = "dotted")
  p=p+geom_hline(yintercept = drug_cutoff,linetype = "dotted")
  p=p+geom_hline(yintercept = -drug_cutoff,linetype = "dotted")
  p=p+geom_abline(slope=slope, intercept=+intercept,linetype = "dotted")
  p=p+geom_abline(slope=slope, intercept=-intercept,linetype = "dotted")
  p=p+labs(x="Control Beta Score",y="Treatment Beta Score",title=main)
  p=p+theme(legend.position="bottom")+theme(legend.title=element_blank())
  p=p+guides(col = guide_legend(ncol = 3,byrow =T))
  p=p+xlim(x_min,x_max)+ylim(y_min,y_max)
  p=p+annotate("text",color="black",
               label=paste("",as.character(dim(beta[beta$group=="Group2",])[1]),sep=""),
               x=0, y=y_max-0.01,vjust=0,hjust = 0.5)
  p=p+annotate("text",color="black",
               label=paste("",as.character(dim(beta[beta$group=="Group4",])[1]),sep=""),
               x=0, y=y_min+0.01,vjust=1,hjust = 0.5)
  p=p+annotate("text",color="black",
               label=paste("",as.character(dim(beta[beta$group=="Group3",])[1]),sep=""),
               x=x_max-0.01, y=0,vjust=0.5,hjust = 1)
  p=p+annotate("text",color="black",
               label=paste("",as.character(dim(beta[beta$group=="Group1",])[1]),sep=""),
               x=x_min+0.01, y=0,vjust=0.5,hjust = 0)
  p=suppressWarnings(ggExtra::ggMarginal(p, type="histogram",bins=50))

  # grid.arrange(p,ncol = 1)
  if(!is.null(filename)){
    write.table(beta, file.path(out.dir,paste0("Square9_detail_",filename,".txt")),
                sep="\t", row.names = F,col.names = T,quote=F)
    ggsave(p,filename=file.path(out.dir,paste0("Square9_scatter_",filename,".png")),
           units="in",width=520/100,height=500/100)
  }
  return(list(p=p,dd1=beta))
}
