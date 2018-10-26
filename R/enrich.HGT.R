#' Do enrichment analysis using Hypergeometric test
#'
#' @docType methods
#' @name enrich.HGT
#' @rdname enrich.HGT
#' @aliases Hypergeometric
#'
#' @param geneList A numeric vector with gene as names.
#' @param universe A character vector, specifying the backgound genelist, default is whole genome.
#' @param keytype "Entrez" or "Symbol".
#' @param type Geneset category for testing, "MsigDB" (default), "KEGG" or "Others" (`gmtpath` is required).
#' @param organism A character, specifying organism, such as "hsa" or "Human"(default), and "mmu" or "Mouse"
#' @param pvalueCutoff Pvalue cutoff.
#' @param pAdjustMethod One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param limit A integer vector of length two, specifying the limit of geneset size.
#' @param gmtpath The path to gmt file.
#'
#' @return A enrichResult instance.
#'
#' @author Wubing Zhang
#'
#' @seealso \code{\link{enrich.GSE}}
#' @seealso \code{\link{enrich.ORT}}
#' @seealso \code{\link{enrichment_analysis}}
#' @seealso \code{\link[DOSE]{enrichResult-class}}
#'
#' @examples
#' data(geneList, package = "DOSE")
#' genes <- geneList[1:100]
#' enrichRes <- enrich.HGT(genes)
#' head(enrichRes@result)
#'
#' @import DOSE
#' @importFrom data.table fread
#'
#' @export

enrich.HGT = function(geneList, universe=NULL, keytype = "Entrez", type="MsigDB", organism='hsa',
                      pvalueCutoff = 0.05, pAdjustMethod = "BH", limit = c(3, 50), gmtpath = NA){
  requireNamespace("clusterProfiler", quietly=TRUE) || stop("need clusterProfiler package")
  requireNamespace("data.table", quietly=TRUE) || stop("need data.table package")
  # download Kegg data
  organism = getOrg(organism)$org
  #============
  if(type=="KEGG"){
    # Prepare pathway annotations
    pathwayFiles <- c(file.path(system.file("extdata", package = "MAGeCKFlute"),
                                paste0("pathways_", organism)),
                      file.path(system.file("extdata", package = "MAGeCKFlute"),
                                paste0("gene2path_", organism)))
    if(!all(file.exists(pathwayFiles))){
      ## Download pathway annotation
      remfname1 <- paste0("http://rest.kegg.jp/link/pathway/", organism)
      remfname2 <- paste0("http://rest.kegg.jp/list/pathway/", organism)
      download.file(remfname1, pathwayFiles[1], quiet = TRUE)
      download.file(remfname2, pathwayFiles[2], quiet = TRUE)
    }
    ## Read and preprocess pathway annotation
    gene2path = data.table::fread(pathwayFiles[1], header = FALSE, showProgress = FALSE)
    names(gene2path)=c("Gene","PathwayID")
    gene2path$PathwayID=gsub("path:","",gene2path$PathwayID)
    gene2path$Gene=gsub(paste0(organism,":"), "", gene2path$Gene)
    pathways = data.table::fread(pathwayFiles[2], header = FALSE, showProgress = FALSE)
    names(pathways) = c("PathwayID","PathwayName")
    pathways$PathwayID = gsub("path:","",pathways$PathwayID)
    pathways = as.data.frame(pathways)
    pathways$PathwayName = gsub(" - .*", "", pathways$PathwayName)
    rownames(pathways) = pathways$PathwayID
    gene2path$PathwayName = pathways[gene2path$PathwayID, "PathwayName"]
  }
  if(type %in% c("MsigDB", "Others")){
    msigdb = file.path(system.file("extdata", package = "MAGeCKFlute"),
                       paste0("MsigDB.all.v6.2.entrez.", organism, ".gmt.gz"))
    if(type == "MsigDB"){ gmtpath = gzfile(msigdb) }
    gene2path = readGMT(gmtpath, limit = limit)
    names(gene2path) = c("Gene","PathwayID", "PathwayName")
    gene2path = gene2path[!is.na(gene2path$Gene), ]
    pathways = data.frame(PathwayID = unique(gene2path$PathwayID),
                          PathwayName = unique(gene2path$PathwayName))
  }
  ##==========
  if(keytype == "Symbol"){
    allsymbol = names(geneList)
    gene = TransGeneID(allsymbol, "Symbol", "Entrez", organism = organism)
    geneList = geneList[!duplicated(gene)]
    names(geneList) = gene[!duplicated(gene)]
    universe = TransGeneID(universe, "Symbol", "Entrez", organism = organism)
  }
  universe = universe[!is.na(universe)]
  gene = names(geneList)
  allsymbol = TransGeneID(gene, "Entrez", "Symbol", organism = organism)
  gene2path$Gene = as.character(gene2path$Gene)
  if(!is.null(universe)){
    universe = unique(as.character(universe))
  }else{
    universe = unique(c(gene, gene2path$Gene))
  }

  # Start to do the hypergeometric test
  m = length(gene)
  n = length(universe) - m
  res = data.frame(ID = c(), Description = c(), GeneRatio = c(), BgRatio = c(),
                   pvalue = c(), geneID = c(),geneName = c(), Count = c())
  kk = 1
  for(c in pathways$PathwayID){
    pathwayGene = gene2path$Gene[gene2path$PathwayID%in%c]
    idx1 = universe %in% pathwayGene
    k = length(universe[idx1])
    idx2 = universe %in% gene
    q = length(universe[idx1&idx2])
    geneID = paste(universe[idx1&idx2], collapse = "/")

    if(k>=limit[1] & k<=limit[2] & q>0){
      pvalue = phyper(q,m,n,k,lower.tail = FALSE)
      res[kk, "ID"] = c
      res[kk, "Description"] = pathways$PathwayName[pathways$PathwayID %in% c]
      res[kk, "GeneRatio"] = paste(q,length(which(idx1)), sep="/")
      res[kk, "BgRatio"] = paste(length(which(idx1)), length(universe), sep="/")
      res[kk, "pvalue"] = pvalue
      np = -log10(pvalue); np[np>10] = 10
      res[kk,"NES"] = sum(geneList[universe[idx1&idx2]]) * np / (sum(idx1&idx2)^0.5)
      res[kk,"geneID"] = geneID
      SYMBOL = allsymbol[universe[idx1&idx2]]
      geneName = paste(SYMBOL, collapse = "/")
      res[kk,"geneName"] = geneName
      res[kk,"Count"] = q
      kk = kk + 1
    }
  }
  res$p.adjust = p.adjust(res$pvalue,pAdjustMethod)
  res$nLogpvalue = -log10(res$p.adjust)
  idx = which(res$pvalue<=pvalueCutoff & res$p.adjust<=pvalueCutoff)
  if(length(idx)>0){
    res = res[idx, ]
    res = res[order(res$NES, decreasing = TRUE), ]
    res = res[, c("ID", "Description", "NES", "pvalue", "p.adjust", "GeneRatio",
                 "BgRatio", "geneID", "geneName", "Count")]
  }else res=data.frame()
  # Create enrichResult object
  new("enrichResult",
      result         = res,
      pvalueCutoff   = pvalueCutoff,
      pAdjustMethod  = pAdjustMethod,
      organism       = organism,
      ontology       = type, ## as.character(x$Category[1]),
      gene           = as.character(gene),
      keytype        = keytype)
}

