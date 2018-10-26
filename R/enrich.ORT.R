#' Do enrichment analysis using over-representation test
#'
#' Do enrichment analysis using over-representation test
#'
#' @docType methods
#' @name enrich.ORT
#' @rdname enrich.ORT
#' @aliases enrichORT
#'
#' @param geneList A numeric vector with gene as names.
#' @param universe A character vector, specifying the backgound genelist, default is whole genome.
#' @param keytype "Entrez" or "Symbol".
#' @param type Geneset category for testing, "MsigDB" (default), "KEGG", "BP", "CC", "MF", "DO",
#' "MKEGG", "NCG" or "Others" (`gmtpath` is required).
#' @param organism A character, specifying organism, such as "hsa" or "Human"(default), and "mmu" or "Mouse"
#' @param pvalueCutoff Pvalue cutoff.
#' @param pAdjustMethod One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param limit A integer vector of length two, specifying the limit of geneset size.
#' @param gmtpath The path to gmt file.
#'
#' @return A enrichedResult instance.
#'
#' @author Wubing Zhang
#'
#' @seealso \code{\link{enrich.HGT}}
#' @seealso \code{\link{enrich.GSE}}
#' @seealso \code{\link{enrichment_analysis}}
#' @seealso \code{\link[DOSE]{enrichResult-class}}
#'
#' @examples
#' data(geneList, package = "DOSE")
#' genes <- geneList[1:100]
#' enrichedRes <- enrich.ORT(genes)
#' head(enrichedRes@result)
#'
#' @import DOSE
#' @import clusterProfiler
#' @export

enrich.ORT <- function(geneList, universe=NULL, keytype = "Entrez", type="MsigDB", organism='hsa',
                       pvalueCutoff = 0.05, pAdjustMethod = "BH", limit = c(3, 50), gmtpath = NA){
  requireNamespace("clusterProfiler", quietly=TRUE) || stop("need clusterProfiler package")
  requireNamespace("data.table", quietly=TRUE) || stop("need data.table package")
  ## Download pathway annotation
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
  ## Gene ID conversion
  if(keytype == "Symbol"){
    allsymbol = names(geneList)
    gene = TransGeneID(allsymbol, "Symbol", "Entrez", organism = organism)
    geneList = geneList[!duplicated(gene)]
    names(geneList) = gene[!duplicated(gene)]
    universe = TransGeneID(universe, "Symbol", "Entrez", organism = organism)
  }

  ## Basic parameters
  universe = universe[!is.na(universe)]
  gene = names(geneList)
  organism = getOrg(organism)$org
  orgdb = getOrg(organism)$pkg
  #=======================
  if(type %in% c("KEGG", "MsigDB")){
    enrichedRes = enricher(gene, pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod,
                  universe = universe, minGSSize = limit[1], maxGSSize = limit[2],
                  TERM2GENE = gene2path[,c("PathwayID","Gene")], TERM2NAME = pathways)
  }
  if(type %in% c("BP", "CC", "MF")){
    enrichedRes = enrichGO(gene=gene,  universe=universe,  ont = type, OrgDb=orgdb,
                           pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff,
                           minGSSize=limit[1], maxGSSize=limit[2])
  }
  if(type == "DO"){
    enrichedRes = enrichDO(gene=gene,  universe=universe, pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff,
                           minGSSize = limit[1], maxGSSize = limit[2])
  }
  if(type == "MKEGG"){
    enrichedRes = enrichMKEGG(gene = gene,  universe=universe, organism = organism,pAdjustMethod=pAdjustMethod,
                              pvalueCutoff = pvalueCutoff, minGSSize = limit[1], maxGSSize = limit[2])
  }
  if(type == "NCG"){
    enrichedRes = enrichNCG(gene=gene, universe = universe, pAdjustMethod = pAdjustMethod,
                            pvalueCutoff = pvalueCutoff, minGSSize = limit[1], maxGSSize = limit[2])
  }
  if(!is.null(enrichedRes)){
    # loginfo("Add symbol to enrichment results ...")
    geneID = strsplit(enrichedRes@result$geneID, split = "/")
    allsymbol = TransGeneID(gene, "Entrez", "Symbol", organism = organism)
    geneName = lapply(geneID, function(gid){
      SYMBOL = allsymbol[gid]
      paste(SYMBOL, collapse = "/")
    })
    enrichedRes@result$geneName = unlist(geneName)
    enrichedRes@result$NES = unlist(apply(enrichedRes@result, 1, function(x){
      enrichGenes = unlist(strsplit(x["geneID"], split = "/"))
      return(sum(geneList[enrichGenes]))
    }))
    enrichedRes@result$NES = enrichedRes@result$NES * (-log10(enrichedRes@result$p.adjust)) / enrichedRes@result$Count^0.5
    enrichedRes@result = enrichedRes@result[order(enrichedRes@result$NES, decreasing = TRUE), ]
    enrichedRes@result = enrichedRes@result[,c("ID", "Description", "NES", "pvalue", "p.adjust",
                                               "GeneRatio", "BgRatio", "geneID", "geneName", "Count")]
  }
  return(enrichedRes)
}

