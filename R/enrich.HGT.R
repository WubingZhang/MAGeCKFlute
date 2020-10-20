#' Do enrichment analysis using hypergeometric test
#'
#' @docType methods
#' @name enrich.HGT
#' @rdname enrich.HGT
#'
#' @param geneList A numeric vector with gene as names
#' @param keytype "Entrez", "Ensembl", or "Symbol"
#' @param type Molecular signatures for testing, available datasets include
#' Pathway (KEGG, REACTOME, C2_CP), GO (GOBP, GOCC, GOMF),
#' MSIGDB (C1, C2 (C2_CP (C2_CP_PID, C2_CP_BIOCARTA), C2_CGP),
#' C3 (C3_MIR, C3_TFT), C4, C6, C7, HALLMARK)
#' and Complex (CORUM). Any combination of them are also accessible
#' (e.g. 'GOBP+GOMF+KEGG+REACTOME')
#' @param organism 'hsa' or 'mmu'
#' @param pvalueCutoff FDR cutoff
#' @param limit A two-length vector, specifying the minimal and
#' maximal size of gene sets for enrichent analysis
#' @param universe A character vector, specifying the backgound genelist, default is whole genome
#' @param gmtpath The path to customized gmt file
#' @param verbose Boolean
#' @param ... Other parameter
#'
#' @return An enrichResult instance.
#'
#' @author Wubing Zhang
#'
#' @seealso \code{\link{enrich.GSE}}
#' @seealso \code{\link{enrich.ORT}}
#' @seealso \code{\link{EnrichAnalyzer}}
#' @seealso \code{\link[DOSE]{enrichResult-class}}
#'
#' @examples
#' data(geneList, package = "DOSE")
#' genes <- geneList[1:300]
#' enrichRes <- enrich.HGT(genes, type = "KEGG", keytype = "entrez")
#' head(slot(enrichRes, "result"))
#'
#' @export

enrich.HGT = function(geneList,
                      keytype = "Symbol",
                      type = "GOBP",
                      organism = 'hsa',
                      pvalueCutoff = 0.25,
                      limit = c(2, 200),
                      universe = NULL,
                      gmtpath = NULL,
                      verbose = TRUE,
                      ...){

  requireNamespace("clusterProfiler", quietly=TRUE) || stop("need clusterProfiler package")

  ## Prepare gene set annotation
  gene2path = gsGetter(gmtpath, type, limit, organism)
  idx = duplicated(gene2path$PathwayID)
  pathways = data.frame(PathwayID = gene2path$PathwayID[!idx],
                        PathwayName = gene2path$PathwayName[!idx],
                        stringsAsFactors = FALSE)

  ## Gene ID conversion
  if(tolower(keytype) != "entrez"){
    allsymbol = names(geneList)
    gene = TransGeneID(allsymbol, keytype, "entrez", organism = organism)
    idx = duplicated(gene)|is.na(gene)
    allsymbol = allsymbol[!idx]; names(allsymbol) = gene[!idx]
    geneList = geneList[!idx]; names(geneList) = gene[!idx]
  }else{
    gene = names(geneList)
    allsymbol = TransGeneID(gene, "Entrez", "Symbol", organism = organism)
  }
  if(!is.null(universe)){
    universe = TransGeneID(universe, keytype, "Entrez", organism = organism)
    universe = universe[!is.na(universe)]
  }else{
    universe = unique(gene2path$Gene)
  }
  gene = names(geneList)

  ## Start to do the hypergeometric test ##
  HGT <- function(pid){
    pGene = gene2path$Gene[gene2path$PathwayID==pid]
    idx1 = universe %in% pGene
    k = length(universe[idx1])
    idx2 = universe %in% gene
    q = length(universe[idx1&idx2])
    geneID = paste(universe[idx1&idx2], collapse = "/")
    retr <- list(ID = NA, Description = NA, NES = NA,
                 pvalue = NA, GeneRatio = NA, BgRatio = NA,
                 geneID = NA, geneName = NA, Count = NA)
    if(k>=limit[1] & k<=limit[2] & q>0){
      pvalue = phyper(q, m, n, k, lower.tail = FALSE)
      retr <- list(ID = pid, Description = pathways$PathwayName[pathways$PathwayID==pid],
                   NES = mean(geneList[universe[idx1&idx2]]) * sum(idx1&idx2)^0.6,
                   pvalue = pvalue, GeneRatio = paste(q, sum(idx2), sep="/"),
                   BgRatio = paste(length(pGene), length(universe), sep="/"),
                   geneID = geneID,
                   geneName = paste(allsymbol[universe[idx1&idx2]],
                                                     collapse = "/"), Count = q)
    }
    return(retr)
  }

  ## Test using above function ##
  len = length(unique(intersect(gene, gene2path$Gene)))
  if(verbose) message("\t", len, " genes are mapped ...")
  m = length(gene)
  n = length(universe) - m
  res = sapply(pathways$PathwayID, HGT)
  res = as.data.frame(t(res), stringsAsFactors = FALSE)
  res = res[!is.na(res$ID), ]
  if(nrow(res)>0){
    res[, c(1:2, 5:8)] = matrix(unlist(res[, c(1:2, 5:8)]), ncol = 6)
    res[, c(3:4, 9)] = matrix(unlist(res[, c(3:4, 9)]), ncol = 3)
    res$p.adjust = p.adjust(res$pvalue, "BH")
    res$nLogpvalue = -log10(res$p.adjust)
    idx = which(res$p.adjust<=pvalueCutoff)
    if(length(idx)>0){
      res = res[idx, ]
      idx = c("ID", "Description", "NES", "pvalue", "p.adjust",
              "GeneRatio", "BgRatio", "geneID", "geneName", "Count")
      res = res[, idx]
    }else res=data.frame()
  }
  res = res[order(res$pvalue), ]
  ## Create enrichResult object ##
  new("enrichResult",
      result         = res,
      pvalueCutoff   = pvalueCutoff,
      pAdjustMethod  = "BH",
      organism       = organism,
      ontology       = type,
      gene           = as.character(gene),
      keytype        = keytype)
}

