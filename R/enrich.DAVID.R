# A job with more than 3000 genes to generate gene or term cluster report will not be handled by DAVID due to resource limit.
# No more than 200 jobs in a day from one user or computer.
# DAVID Team reserves right to suspend any improper uses of the web service without notice.
#' Do enrichment analysis using DAVID
#'
#' an update version of DAVIDWebService to do enrichment analysis
#'
#' @docType methods
#' @name enrich.DAVID
#' @rdname enrich.DAVID
#' @aliases enrichDAVID
#'
#' @param gene Character vector, specifying the genelist to do enrichment analysis.
#' @param universe Character vector, specifying the backgound genelist, default is whole genome.
#' @param david.user Character, specifying a valid DAVID user account.
#' @param idType Character, indicating the gene id type of input genelist, such as "ENTREZ_GENE_ID"(default).
#' @param minGSSize Minimal size of each geneSet for testing.
#' @param maxGSSize Maximal size of each geneSet for analyzing.
#' @param annotation Geneset category for testing, GOTERM_BP_FAT(default).
#' @param pvalueCutoff Pvalue cutoff.
#' @param pAdjustMethod One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param qvalueCutoff Qvalue cutoff.
#'
#' @return A enrichResult instance.
#'
#' @author Wubing Zhang
#'
#' @note This function depends on network and DAVID account, so don't show in the vignette.
#' The source can be found by typing \code{MAGeCKFlute:::enrich.DAVID}
#' or \code{getMethod("enrich.DAVID")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/enrich.DAVID.R}
#' Users should find it easy to customize this function.
#'
#' @seealso \code{\link{enrich.HGT}}
#' @seealso \code{\link{enrich.GOstats}}
#' @seealso \code{\link{enrich.GSE}}
#' @seealso \code{\link{enrich.ORT}}
#' @seealso \code{\link{enrichment_analysis}}
#' @seealso \code{\link[DOSE]{enrichResult-class}}
#'
#' @examples
#' \dontrun{
#'  data(MLE_Data)
#'  universe = TransGeneID(MLE_Data$Gene, "SYMBOL", "ENTREZID", organism = "hsa")
#'  genes = universe[1:50]
#' 	# Before running this example, you need to have a david account.
#' 	enrichRes <- enrich.DAVID(genes, universe, david.user="david.user@edu.com")
#' 	head(enrichRes@result)
#' }
#'
#'
#' @export


enrich.DAVID <- function(gene, universe=NULL, david.user, idType="ENTREZ_GENE_ID",
                         minGSSize = 2, maxGSSize = 500, annotation  = "GOTERM_BP_FAT",
                         pvalueCutoff  = 1, pAdjustMethod = "BH", qvalueCutoff= 0.2){

  ## user:ma.tongji@gmail.com
  david.pkg <- "RDAVIDWebService"
  pkgs <- installed.packages()[,1]
  if (! david.pkg %in% pkgs) {
    stop("Install RDAVIDWebService package before using enrichDAVID...")
  }
  # requireNamespace(david.pkg)
  Count <- List.Total <- Pop.Hits <- Pop.Total <- NULL

  pAdjustMethod <- match.arg(pAdjustMethod, c("bonferroni", "BH"))

  DAVIDWebService <- eval(parse(text="DAVIDWebService"))
  addList <- eval(parse(text="addList"))
  setAnnotationCategories <- eval(parse(text="setAnnotationCategories"))
  getFunctionalAnnotationChart <- eval(parse(text="getFunctionalAnnotationChart"))
  getSpecieNames <- eval(parse(text="getSpecieNames"))
  getIdTypes <- eval(parse(text="getIdTypes"))

  david <- DAVIDWebService$new(email=david.user,
                               url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

  ## addList will throw error if idType is not match.
  ## use match.arg to check before addList make it more readable

  idType <- match.arg(idType, getIdTypes(david))

  david.res <- addList(david, gene, idType=idType, listName="clusterProfiler", listType="Gene")
  if(!is.null(universe)){
    david.res <- addList(david, universe, idType=idType, listName="User_assign", listType="Background")
  }

  if (david.res$inDavid == 0) {
    stop("All id can not be mapped. Please check 'idType' parameter...")
  }

  setAnnotationCategories(david, annotation)
  x <- getFunctionalAnnotationChart(david, threshold=1, count=minGSSize)

  if (length(x@.Data) == 0) {
    warning("No significant enrichment found...")
    return(NULL)
  }

  term <- x$Term
  if (length(grep("~", term[1])) == 0) {
    sep <- ":"
  } else {
    sep <- "~"
  }
  term.list <- sapply(term, function(y) strsplit(y, split=sep))
  term.df <- do.call("rbind", term.list)
  ID <- term.df[,1]
  Description <- term.df[,2]
  GeneRatio <- with(x, paste(Count, List.Total, sep="/"))
  BgRatio <- with(x, paste(Pop.Hits, Pop.Total, sep="/"))
  Over <- data.frame(ID          = ID,
                     Description = Description,
                     GeneRatio   = GeneRatio,
                     BgRatio     = BgRatio,
                     pvalue      = x$PValue,
                     stringsAsFactors = FALSE)
  row.names(Over) <- ID

  if (pAdjustMethod == "bonferroni") {
    Over$p.adjust <- x$Bonferroni
  } else {
    Over$p.adjust <- x$Benjamini
  }

  qobj <- tryCatch(qvalue(p=Over$pvalue, lambda=0.05, pi0.method="bootstrap"),
                   error=function(e) NULL)
  if (class(qobj) == "qvalue") {
    qvalues <- qobj$qvalues
  } else {
    qvalues <- NA
  }
  Over$qvalue <- qvalues
  Over$geneID <- gsub(",\\s*", "/", x$Genes)
  gc <- unlist(strsplit(Over$geneID, "/"))

  org <- getSpecieNames(david)
  org <- gsub("\\(.*\\)", "", org)

  if(org=="Homo sapiens"){
    SYMBOL = suppressMessages(eg2id(gc, "SYMBOL", pkg.name = getOrg("hsa")$lib)[, "SYMBOL"])
  }else{
    SYMBOL = suppressMessages(eg2id(gc, "SYMBOL", pkg.name = getOrg("mmu")$lib)[, "SYMBOL"])
  }
  geneName=paste(SYMBOL, collapse = "/")
  Over$geneName <- geneName

  Over$Count <- x$Count

  Over <- Over[ Over$pvalue <= pvalueCutoff, ]
  Over <- Over[ Over$p.adjust <= pvalueCutoff, ]
  if (! any(is.na(Over$qvalue))) {
    Over <- Over[ Over$qvalue <= qvalueCutoff, ]
  }




  if (!is.na(maxGSSize) || !is.null(maxGSSize)) {
    idx <- as.numeric(sub("/\\d+", "", Over$BgRatio)) <= maxGSSize
    Over <- Over[idx,]
  }

  new("enrichResult",
      result         = Over,
      pvalueCutoff   = pvalueCutoff,
      pAdjustMethod  = pAdjustMethod,
      organism       = org,
      ontology       = annotation, ## as.character(x$Category[1]),
      gene           = as.character(gene),
      keytype        = idType)
}

