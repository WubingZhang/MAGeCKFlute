enrichment_plot=function(enrichment,plotTitle=NULL,gridColour="blue",termNum=10,charLength=50){

  if(is.null(enrichment) || nrow(enrichment)==0){
    p1=ggplot()
    p1=p1+geom_text(aes(x=0,y=0,label="No enriched terms"),size=6)
    p1=p1+labs(title=plotTitle)
    p1=p1+theme(plot.title = element_text(size=10))
    p1=p1+theme_void()
    p1=p1+theme(plot.title = element_text(hjust = 0.5))
    return(p1)
  }

  if(nrow(enrichment)>=termNum){
    enrichment=enrichment[1:termNum,]
  }

  #The column of Description, ID, p.adjust, and Count are neccessary.
  enrichment$logP = -log10(enrichment$p.adjust)
  enrichment = enrichment[!is.na(enrichment$ID),]
  enrichment=enrichment[!duplicated(enrichment$Description),]
  enrichment = enrichment[order(enrichment$logP,decreasing = TRUE), ]

  #normalize term description
  {
    terms=as.character(enrichment$Description)
    terms=lapply(terms,function(x,k){
      x=as.character(x)
      x=paste(toupper(substr(x,1,1)),substr(x,2,nchar(x)),sep="")
      if(nchar(x)>k){
        x=substr(x,start=1,stop=k)
        x=gsub("(\\w+)$", "...", x)
      }
      return(x)
    },charLength)
    enrichment$Description=do.call(rbind,terms)
  }
  idx = duplicated(enrichment$Description)
  enrichment=enrichment[!idx,]
  enrichment$Name=factor(enrichment$Description,levels=rev(enrichment$Description))

  p1=ggplot(data=enrichment,aes(x=logP,y=Name,size=Count,colour = factor(Count)))
  p1=p1+geom_point()
  p1=p1+guides(color= FALSE)
  p1 <- p1+theme(panel.grid.major=element_line(colour=gridColour),
                 panel.grid.minor=element_blank(),
                 panel.background=element_blank())
  p1 <- p1+xlab("-log10(Adjust.pvalue)")+ylab("")
  p1 <- p1+labs(title=plotTitle)
  p1 <- p1+theme(axis.text.x=element_text(size=6, face="plain", colour='black'))
  p1 <- p1+theme(axis.text.y=element_text(size=6, face="plain", colour='black'))
  p1=p1+theme(legend.position="bottom")+theme(plot.title = element_text(hjust = 0.5,size=10))
  #p1
  return(p1)
}

