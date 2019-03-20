#' Do enrichment analysis using Hypergeometric test
#'
#' @docType methods
#' @name enrich.HGT
#' @rdname enrich.HGT
#' @aliases Hypergeometric
#'
#' @param geneList A numeric vector with gene as names.
#' @param keytype "Entrez" or "Symbol".
#' @param type Geneset category for testing, one of 'CORUM', 'CPX' (ComplexPortal),
#' 'GOBP', 'GOMF', 'GOCC', 'KEGG', 'BIOCARTA', 'REACTOME', 'WikiPathways', 'EHMN', 'PID',
#' or any combination of them (e.g. 'GOBP+GOMF+CORUM'), or 'All' (all categories).
#' @param organism 'hsa' or 'mmu'.
#' @param pvalueCutoff Pvalue cutoff.
#' @param limit A two-length vector (default: c(3, 50)), specifying the minimal and
#' maximal size of gene sets for enrichent analysis.
#' @param universe A character vector, specifying the backgound genelist, default is whole genome.
#' @param gmtpath The path to customized gmt file.
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
#' genes <- geneList[1:300]
#' enrichRes <- enrich.HGT(genes, type = "KEGG")
#' head(slot(enrichRes, "result"))
#'
#' @import DOSE
#' @importFrom data.table fread
#'
#' @export

enrich.HGT = function(geneList, keytype = "Entrez",
                      type = "CORUM+GOBP+GOMF+RECTOME+KEGG",
                      organism = 'hsa', pvalueCutoff = 0.05,
                      limit = c(3, 50), universe = NULL, gmtpath = NA){
  requireNamespace("clusterProfiler", quietly=TRUE) || stop("need clusterProfiler package")
  requireNamespace("data.table", quietly=TRUE) || stop("need data.table package")

  ## Prepare gene set annotation
  if(is.na(gmtpath)){
    msigdb = file.path(system.file("extdata", package = "MAGeCKFlute"),
                       paste0(organism, "_msig_entrez.gmt.gz"))
    gmtpath = gzfile(msigdb)
  }
  gene2path = ReadGMT(gmtpath, limit = limit)
  close(gmtpath)
  names(gene2path) = c("Gene","PathwayID", "PathwayName")
  gene2path$PathwayName = paste0(toupper(substr(gene2path$PathwayName, 0, 1)),
                                 tolower(substr(gene2path$PathwayName, 2,
                                                nchar(gene2path$PathwayName))))
  if(type != "All"){
    type = unlist(strsplit(type, "\\+"))
    idx = toupper(gsub("_.*", "", gene2path$PathwayID)) %in% toupper(type)
    gene2path = gene2path[idx, ]
  }
  gene2path = gene2path[!is.na(gene2path$Gene), ]
  idx = duplicated(gene2path$PathwayID)
  pathways = data.frame(PathwayID = gene2path$PathwayID[!idx],
                        PathwayName = gene2path$PathwayName[!idx],
                        stringsAsFactors = FALSE)

  ## Gene ID conversion
  if(keytype != "Entrez"){
    allsymbol = names(geneList)
    gene = TransGeneID(allsymbol, keytype, "Entrez", organism = organism)
    geneList = geneList[!duplicated(gene)]
    names(geneList) = gene[!duplicated(gene)]
    universe = TransGeneID(universe, keytype, "Entrez", organism = organism)
    universe = universe[!is.na(universe)]
  }
  gene = names(geneList)
  allsymbol = TransGeneID(gene, "Entrez", "Symbol", organism = organism)
  gene2path$Gene = as.character(gene2path$Gene)
  if(!is.null(universe)){
    universe = unique(as.character(universe))
  }else{
    universe = unique(c(gene, gene2path$Gene))
  }

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
    if(k>=limit[1] & k<=limit[2] & q>2){
      pvalue = phyper(q, m, n, k, lower.tail = FALSE)
      retr <- list(ID = pid, Description = pathways$PathwayName[pathways$PathwayID==pid],
                   NES = mean(geneList[universe[idx1&idx2]]) * sum(idx1&idx2)^0.6,
                   pvalue = pvalue, GeneRatio = paste(q, sum(idx1), sep="/"),
                   BgRatio = paste(sum(idx1), length(pGene), sep="/"),
                   geneID = geneID, geneName = paste(allsymbol[universe[idx1&idx2]],
                                                     collapse = "/"), Count = q)
    }
    return(retr)
  }

  ## Test using above function ##
  len = length(unique(intersect(gene, gene2path$Gene)))
  message("\t", len, " genes are mapped ...")
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

