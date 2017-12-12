#' Do enrichment analysis using Hypergeometric test
#'
#' Do enrichment analysis using Hypergeometric test
#'
#' @docType methods
#' @name enrich.HGT
#' @rdname enrich.HGT
#' @aliases Hypergeometric
#'
#' @param gene a character vector, specifying the genelist to do enrichment analysis.
#' @param universe a character vector, specifying the backgound genelist, default is whole genome.
#' @param type geneset category for testing, KEGG(default).
#' @param organism a character, specifying organism, such as "hsa" or "Human"(default), and "mmu" or "Mouse"
#' @param pvalueCutoff pvalue cutoff.
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param minGSSize minimal size of each geneSet for testing.
#' @param maxGSSize maximal size of each geneSet for analyzing.
#'
#' @return A enrichResult instance.
#'
#' @author Feizhen Wu
#'
#' @note  See the vignette for an example of enrichment analysis using hypergemetric test
#' The source can be found by typing \code{MAGeCKFlute:::enrich.HGT}
#' or \code{getMethod("enrich.HGT")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/enrich.HGT.R}
#' Users should find it easy to customize this function.
#'
#' @seealso \code{\link{enrich.GOstats}}
#' @seealso \code{\link{enrich.DAVID}}
#' @seealso \code{\link{enrich.GSE}}
#' @seealso \code{\link{enrich.ORT}}
#' @seealso \code{\link{enrichment_analysis}}
#' @seealso \code{\link[DOSE]{enrichResult-class}}
#'
#' @examples
#' data(MLE_Data)
#' universe = TransGeneID(MLE_Data$Gene, "SYMBOL", "ENTREZID", organism = "hsa")
#' genes = TransGeneID(Core_Essential[1:200], "SYMBOL", "ENTREZID", organism = "hsa")
#' enrichRes <- enrich.HGT(genes, universe)
#' head(enrichRes@result)
#'
#' @importFrom data.table fread
#' @importFrom pathological temp_dir
#'
#' @export


enrich.HGT = function(gene, universe, type="KEGG", organism='hsa', pvalueCutoff = 1,
                      pAdjustMethod = "BH", minGSSize = 2, maxGSSize = 500){
  requireNamespace("data.table", quietly=TRUE) || stop("need data.table package")
  requireNamespace("pathological", quietly=TRUE) || stop("need pathological package")

  gene = unique(as.character(gene))
  universe = unique(as.character(universe))
  # download Kegg data
  organism = getOrg(organism)$org
  pathwayFiles <- c(file.path(temp_dir(), paste0("pathways_", organism)),
                    file.path(temp_dir(), paste0("gene2path_", organism)))

  if(!all(file.exists(pathwayFiles))){
    gene2path=fread(paste0("http://rest.kegg.jp/link/pathway/",organism),
                    header = FALSE, showProgress = FALSE)
    names(gene2path)=c("EntrezID","PathwayID")
    gene2path$PathwayID=gsub("path:","",gene2path$PathwayID)
    gene2path$EntrezID=gsub(paste0(organism,":"),"",gene2path$EntrezID)

    pathways=fread(paste0("http://rest.kegg.jp/list/pathway/",organism),
                   header = FALSE, showProgress = FALSE)
    names(pathways)=c("PathwayID","PathwayName")

    pathways$PathwayID=gsub("path:","",pathways$PathwayID)
    pathways$PathwayName=gsub(" - .*", "", pathways$PathwayName)

    write.table(pathways, pathwayFiles[1], sep="\t", row.names = FALSE)
    write.table(gene2path, pathwayFiles[2], sep="\t", row.names = FALSE)
  }else{
    pathways=read.table(pathwayFiles[1], sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    gene2path=read.table(pathwayFiles[2], sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  }

  #============
  loginfo(paste('Running KEGG patwhay for list of entrezIDs'))
  if(type=="KEGG"){
    m=length(gene)
    n=length(universe)-m
    res=data.frame(ID=c(),Description=c(),GeneRatio=c(),BgRatio=c(),pvalue=c(),
                   geneID=c(),geneName=c(), Count=c())
    kk=1
    for(c in pathways$PathwayID){
      pathwayEntrezID=gene2path$EntrezID[gene2path$PathwayID%in%c]
      idx1=universe %in% pathwayEntrezID
      k=length(universe[idx1])
      idx2=universe %in% gene
      q=length(universe[idx1&idx2])
      geneID=paste(universe[idx1&idx2],collapse = "/")

      if(k>=minGSSize&k<=maxGSSize&q>0){
        pvalue=phyper(q,m,n,k,lower.tail = FALSE)
        res[kk,"ID"]=c
        res[kk,"Description"]=pathways$PathwayName[pathways$PathwayID %in% c]
        res[kk,"GeneRatio"]=paste(q,length(which(idx1)),sep="/")
        res[kk,"BgRatio"]=paste(length(which(idx1)),length(gene2path$EntrezID),sep="/")
        res[kk,"pvalue"]=pvalue
        res[kk,"geneID"]=geneID

        SYMBOL = TransGeneID(universe[idx1&idx2], "ENTREZID", "SYMBOL", organism = organism)
        geneName = paste(SYMBOL, collapse = "/")
        res[kk,"geneName"]=geneName
        res[kk,"Count"]=q

        kk=kk+1
      }
    }
    res$p.adjust=p.adjust(res$pvalue,pAdjustMethod)
    res$nLogpvalue = -log10(res$p.adjust)
    idx = which(res$pvalue<=pvalueCutoff & res$p.adjust<=pvalueCutoff)

    if(length(idx)>0){
      res = res[idx, ]
      res = res[order(res$p.adjust),]
      res = res[,c("ID", "Description", "pvalue", "p.adjust", "nLogpvalue", "GeneRatio",
                   "BgRatio", "geneID", "geneName", "Count")]
    }else{
      res=data.frame()
    }

  }
  new("enrichResult",
      result         = res,
      pvalueCutoff   = pvalueCutoff,
      pAdjustMethod  = pAdjustMethod,
      organism       = organism,
      ontology       = type, ## as.character(x$Category[1]),
      gene           = as.character(gene),
      keytype        = "ENTREZID")
}

