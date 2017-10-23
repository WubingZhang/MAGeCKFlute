getOrg <- function(organism, onlyLib=F){
  data(bods)
  res=list()
  ##======
  # Get the mapping from organism to package
  orgMap = as.data.frame(bods[,1:3], stringsAsFactors=F)
  orgMap$species = tolower(orgMap$species)
  orgMap$`kegg code` = tolower(orgMap$`kegg code`)
  orgMap = melt(orgMap, id="package")
  org.package = orgMap$package
  names(org.package) = orgMap$value

  # Get the mapping from organism to alias2eg
  tmp = t(as.data.frame(strsplit(org.package, "\\.")))
  tmp[,3] = rep("egALIAS2EG",nrow(tmp))
  tmp = tmp[,-4]
  org.alias = apply(tmp, 1, paste, collapse=".")

  #==============
  # Get corresponding package and alias2eg dataframe
  organism = tolower(organism)
  if(organism %in% names(org.package)){
    pkg = org.package[organism]
    alias = org.alias[organism]
    idx = which(org.package==pkg)[2]
    res$pkg = pkg   # return pkg
    res$organism = names(org.package)[idx]
  }else{
    stop("No organism found !")
  }

  # Install package
  if(!pkg%in%installed.packages()){
    source("http://bioconductor.org/biocLite.R")
    biocLite(pkg)
  }
  suppressPackageStartupMessages(library(pkg,character.only=T))

  #===============
  # Get mapping data
  if(!onlyLib){
    aliasMapping <- as.data.frame(toTable(get(alias)))
    tmp=suppressWarnings(suppressMessages(as.data.frame(
      bitr(aliasMapping$gene_id, fromType="ENTREZID", toType="SYMBOL", OrgDb=pkg))))
    aliasMapping = merge(aliasMapping,tmp,by.x="gene_id",by.y="ENTREZID")
    aliasMapping = data.frame(gene_id=aliasMapping$gene_id,
                              SYMBOL=c(aliasMapping$SYMBOL, aliasMapping$alias_symbol))
    tmp=suppressWarnings(suppressMessages(as.data.frame(
      bitr(aliasMapping$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb=pkg))))
    tmp1 = merge(aliasMapping,tmp,by.x="SYMBOL",by.y="SYMBOL")
    mapping = tmp1[,c("SYMBOL","ENTREZID")]
    idx1 = duplicated(mapping$SYMBOL)
    idx2 = duplicated(mapping$ENTREZID)
    mapping = mapping[!(idx1&idx2),]
    mapping$SYMBOL = toupper(mapping$SYMBOL)
    res$Symbol_Entrez = mapping
  }

  # return
  return(res)

}
