
pathPng2Pdf <- function(genelist, pathways=c(), organism='hsa', view_allpath=F, title="Group A", sub="Negative control normalized",
                        output=".", path.archive=".", kegg.native = T){
  #====No pathways supplied======================
  if(length(pathways)<1){
    p=ggplot()
    p=p+geom_text(aes(x=0,y=0,label="No enriched terms"),size=6)
    p=p+theme_void()
    print(p)
    return(0)
  }

  if(length(pathways)<4 || view_allpath){
    keggID=pathways
  }else{
    keggID=pathways[1:4]
  }

  loginfo(paste('Starting plot kegg pathways for',sub, title))
  idx=duplicated(genelist$ENTREZID)
  genelist=genelist[!idx,]
  rownames(genelist)=genelist$ENTREZID


  p1 <- KeggPathwayView(gene.data  = genelist[,c("Control","Treatment")], pathway.id = keggID,
                        species=organism, kegg.dir = path.archive, kegg.native = kegg.native)

  #Maybe there are not multi file, but only keggID.pathview.png
  allpngnames=paste0(keggID, ".pathview.multi.png")
  idx = file.exists(allpngnames)
  allpngnames=allpngnames[idx]

  if(length(allpngnames)>0){
    toFile=paste0(output,"/",title,"_",sub,"_",allpngnames)
    boo=file.rename(from=allpngnames,to=toFile)
  }else{boo=FALSE}
  originPng=paste0(keggID, ".png")
  originXML=paste0(keggID, ".xml")
  failMulti=paste0(keggID, ".pathview.png")
  suppressWarnings(file.remove(originPng))
  suppressWarnings(file.remove(originXML))
  suppressWarnings(file.remove(failMulti))

  if(all(boo)){
    pngnames = paste0(output,"/",title,"_",sub,"_",allpngnames)
    idx=file.exists(pngnames)
    pngnames = pngnames[idx]
  }else{pngnames=c()}

  if(length(pngnames)>0){
    thePlots <- lapply (pngnames, function(figure) {
      rasterGrob(readPNG(figure, native = FALSE),interpolate = FALSE)})
  }else
    thePlots = list()

  if(length(thePlots)<4){
    for(i in (length(thePlots)+1):4){
      p1=ggplot()
      p1=p1+geom_text(aes(x=0,y=0,label="No multi pathview figures"),size=6)
      p1=p1+theme_void()
      thePlots[[i]]=p1
    }
  }
  do.call(grid.arrange, c(thePlots[1:2], ncol = 2,top=title,bottom=sub))
  do.call(grid.arrange, c(thePlots[3:4], ncol = 2,top=title,bottom=sub))
  # grid.arrange(thePlots, ncol = 2, top=title0)
}
