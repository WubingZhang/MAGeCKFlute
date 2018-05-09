#' Determine the gene annotation package.
#'
#' Determine the gene annotation package. for specific organism
#'
#' @docType methods
#' @name getOrg
#' @rdname getOrg
#'
#' @param organism Character, KEGG species code, or the common species name, used to determine
#' the gene annotation package. For all potential values check: data(bods); bods. Default org="hsa",
#' and can also be "human" (case insensitive).
#'
#' @return A list containing three elements:
#' \item{organism}{species}
#' \code{pkg}{annotation package name}
#' \code{Symbol_Entrez}{a data frame, mapping between gene symbol and entrez id}
#'
#' @author Wubing Zhang
#'
#' @note
#' The source can be found by typing \code{MAGeCKFlute:::getOrg}
#' or \code{getMethod("getOrg")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/getOrg.R}
#' Users should find it easy to customize this function.
#'
#' @examples
#' ann_pkg = getOrg("human")$pkg
#' print(ann_pkg)
#'
#' @importFrom reshape melt

#'
#' @export

getOrg <- function(organism){
  requireNamespace("reshape", quietly=TRUE) || stop("need reshape package")
  requireNamespace("clusterProfiler", quietly=TRUE) || stop("need clusterProfiler package")
  data(bods)
  res=list()
  ##======
  # Get the mapping from organism to package
  ridx=grep(tolower(organism), tolower(bods[,2:3])) %% nrow(bods)
  if(length(ridx)==0) stop("Wrong organism name!")
  res$org = bods[ridx,3]
  res$pkg = bods[ridx,1]

  tmpFile = file.path(temp_dir(), paste0("map_symbol_entrez_", res$org))
  if(!file.exists(tmpFile)){
    gzfile = c("Homo_sapiens.gene_info.gz", "Bos_taurus.gene_info.gz",
               "Canis_familiaris.gene_info.gz", "Mus_musculus.gene_info.gz",
               "Pan_troglodytes.gene_info.gz", "Rattus_norvegicus.gene_info.gz",
               "Sus_scrofa.gene_info.gz")
    names(gzfile) = c("hsa", "bta", "cfa", "mmu", "ptr", "rno", "ssc")
    locfname <- file.path(temp_dir(), gzfile[res$org])
    if(!file.exists(locfname)){
      remfname <- paste0("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/", gzfile[res$org])
      download.file(remfname, locfname, quiet = TRUE)
    }
    data = read.table(gzfile(locfname), sep = "\t", header = TRUE,
                      quote = "", stringsAsFactors = FALSE, comment.char = "")
    data = data[, c("GeneID", "Symbol", "Synonyms", "dbXrefs", "type_of_gene", "description")]

    data2 = apply(data, 1, function(x){
      Ensembl = "-"; HGNC = "-"
      if(grepl("HGNC", x[4])){
        HGNC = gsub(".*HGNC:", "", x[4])
        HGNC = gsub("\\|.*", "", HGNC)
      }
      if(grepl("Ensembl", x[4])){
        Ensembl = gsub(".*Ensembl:", "", x[4])
        Ensembl = gsub("\\|.*", "", Ensembl)
      }
      tmp = unlist(strsplit(x[3], "[|]"))
      as.vector(rbind(as.integer(x[1]), x[2], tmp, HGNC, Ensembl, x[5], x[6]))
    })
    data2 = matrix(unlist(data2), ncol=7, byrow = TRUE)
    colnames(data2) = c("ENTREZID", "SYMBOL", "Synonyms", "HGNC", "ENSEMBL", "TYPE", "FULLNAME")
    tmp = reshape::melt(as.data.frame(data2, stringsAsFactors = FALSE), measure.vars=c("SYMBOL", "Synonyms"))
    mapping = tmp[, c("ENTREZID", "value", "ENSEMBL", "HGNC", "TYPE", "FULLNAME")]
    colnames(mapping)[2] = "SYMBOL"
    idx = duplicated(paste(mapping$ENTREZID, mapping$SYMBOL, mapping$ENSEMBL,
                           mapping$HGNC, mapping$TYPE, mapping$FULLNAME, sep = "_"))
    mapping = mapping[!idx, ]
    idx = mapping$SYMBOL=="-" & mapping$ENSEMBL=="-" & mapping$HGNC=="-" &
      mapping$TYPE=="-" & mapping$FULLNAME=="-"
    mapping = mapping[!idx, ]
    mapping$SYMBOL = toupper(mapping$SYMBOL)
    res$Symbol_Entrez = mapping
    write.table(mapping, tmpFile, sep = "\t", row.names = FALSE)
  }else{
    res$Symbol_Entrez = read.table(tmpFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  }

  # # Get the mapping from organism to alias2eg
  # tmp = strsplit(res$pkg, "\\.")[[1]]
  # tmp[3] = "egALIAS2EG"
  # tmp = tmp[-4]
  # alias = paste(tmp, collapse = ".")
  # #==============
  # # Get mapping
  # require(res$pkg, character.only = TRUE)
  # tmpFile = file.path(temp_dir(), paste0("map_symbol_entrez_", res$org))
  # if(!file.exists(tmpFile)){
  #   aliasMapping <- as.data.frame(toTable(get(alias)))
  #   tmp=suppressWarnings(suppressMessages(as.data.frame(
  #     bitr(aliasMapping$gene_id, fromType="ENTREZID", toType="SYMBOL", OrgDb=res$pkg))))
  #   aliasMapping = merge(aliasMapping,tmp,by.x="gene_id",by.y="ENTREZID")
  #   aliasMapping = data.frame(gene_id=aliasMapping$gene_id,
  #                             SYMBOL=c(aliasMapping$SYMBOL, aliasMapping$alias_symbol))
  #   tmp=suppressWarnings(suppressMessages(as.data.frame(
  #     bitr(aliasMapping$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb=res$pkg))))
  #   tmp1 = merge(aliasMapping,tmp,by.x="SYMBOL",by.y="SYMBOL")
  #   mapping = tmp1[,c("SYMBOL","ENTREZID")]
  #   idx1 = duplicated(mapping$SYMBOL)
  #   idx2 = duplicated(mapping$ENTREZID)
  #   mapping = mapping[!(idx1&idx2),]
  #   mapping$SYMBOL = toupper(mapping$SYMBOL)
  #   mapping$ENTREZID = as.character(mapping$ENTREZID)
  #   res$Symbol_Entrez = mapping
  #   write.table(mapping, tmpFile, sep = "\t", row.names = FALSE)
  # }else{
  #   res$Symbol_Entrez = read.table(tmpFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # }
  return(res)
}
