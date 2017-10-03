enrich.DAVID <- function(gene, universe=NULL, david.user, idType="ENTREZ_GENE_ID",
                         minGSSize = 2, maxGSSize = 500, annotation  = "GOTERM_BP_FAT",
                         pvalueCutoff  = 0.05, pAdjustMethod = "BH", qvalueCutoff= 0.2){

  # A job with more than 3000 genes to generate gene or term cluster report will not be handled by DAVID due to resource limit.
  # No more than 200 jobs in a day from one user or computer.
  # DAVID Team reserves right to suspend any improper uses of the web service without notice.
  ##' enrichment analysis by DAVID

  loginfo(paste('Running DAVID for list of entrezIDs'))

  Count <- List.Total <- Pop.Hits <- Pop.Total <- NULL

  pAdjustMethod <- match.arg(pAdjustMethod, c("bonferroni", "BH"))

  david.pkg <- "RDAVIDWebService"
  pkgs <- installed.packages()[,1]
  if (! david.pkg %in% pkgs) {
    stop("You should have RDAVIDWebService package installed before using enrichDAVID...")
  }

  require(david.pkg, character.only=TRUE)
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
    SYMBOL = TransGeneID(gc, "ENTREZID", "SYMBOL", "hsa")
  }else{
    SYMBOL = TransGeneID(gc, "ENTREZID", "SYMBOL", "mmu")
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

