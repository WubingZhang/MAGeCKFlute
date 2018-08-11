#' Do enrichment analysis using Hypergeometric test
#'
#' Do enrichment analysis using Hypergeometric test
#'
#' @docType methods
#' @name enrich.HGT
#' @rdname enrich.HGT
#' @aliases Hypergeometric
#'
#' @param gene A character vector, specifying the genelist to do enrichment analysis.
#' @param universe A character vector, specifying the backgound genelist, default is whole genome.
#' @param type Geneset category for testing, KEGG(default).
#' @param organism A character, specifying organism, such as "hsa" or "Human"(default), and "mmu" or "Mouse"
#' @param pvalueCutoff Pvalue cutoff.
#' @param pAdjustMethod One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param minGSSize Minimal size of each geneSet for testing.
#' @param maxGSSize Maximal size of each geneSet for analyzing.
#'
#' @return A enrichResult instance.
#'
#' @author Feizhen Wu
#'
#' @seealso \code{\link{enrich.GOstats}}
#' @seealso \code{\link{enrich.DAVID}}
#' @seealso \code{\link{enrich.GSE}}
#' @seealso \code{\link{enrich.ORT}}
#' @seealso \code{\link{enrichment_analysis}}
#' @seealso \code{\link[DOSE]{enrichResult-class}}
#'
#' @examples
#' data(geneList, package = "DOSE")
#' genes <- names(geneList)[1:100]
#' enrichRes <- enrich.HGT(genes)
#' head(enrichRes@result)
#'
#' @import DOSE
#' @importFrom data.table fread
#'
#' @export


enrich.HGT = function(gene, universe=NULL, type="KEGG", organism='hsa', pvalueCutoff = 0.25,
                      pAdjustMethod = "BH", minGSSize = 2, maxGSSize = 500){
  requireNamespace("data.table", quietly=TRUE) || stop("need data.table package")
  # download Kegg data
  organism = getOrg(organism)$org
  pathwayFiles <- c(file.path(system.file("extdata", package = "MAGeCKFlute"),
                              paste0("pathways_", organism)),
                    file.path(system.file("extdata", package = "MAGeCKFlute"),
                              paste0("gene2path_", organism)))
  if(!all(file.exists(pathwayFiles))){
    ## Download pathway annotation
    remfname1 <- paste0("http://rest.kegg.jp/link/pathway/",organism)
    remfname2 <- paste0("http://rest.kegg.jp/list/pathway/",organism)
    download.file(remfname1, pathwayFiles[1], quiet = TRUE)
    download.file(remfname2, pathwayFiles[2], quiet = TRUE)
  }
  ## Read and preprocess pathway annotation
  gene2path = data.table::fread(pathwayFiles[1], header = FALSE, showProgress = FALSE)
  names(gene2path)=c("EntrezID","PathwayID")
  gene2path$PathwayID=gsub("path:","",gene2path$PathwayID)
  gene2path$EntrezID=gsub(paste0(organism,":"),"",gene2path$EntrezID)
  pathways=data.table::fread(pathwayFiles[2], header = FALSE, showProgress = FALSE)
  names(pathways)=c("PathwayID","PathwayName")
  pathways$PathwayID=gsub("path:","",pathways$PathwayID)
  pathways$PathwayName=gsub(" - .*", "", pathways$PathwayName)
  ##==========
  gene = unique(as.character(gene))
  allsymbol = TransGeneID(gene, "Entrez", "Symbol", organism = organism)
  if(!is.null(universe)){universe = unique(as.character(universe))
  }else{universe = unique(c(gene, gene2path$EntrezID))}

  #============
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

        SYMBOL = allsymbol[universe[idx1&idx2]]
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

