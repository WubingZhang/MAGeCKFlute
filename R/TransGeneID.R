#' Gene ID conversion between ENTREZID and SYMBOL
#'
#' @docType methods
#' @name TransGeneID
#' @rdname TransGeneID
#' @aliases transGeneID
#'
#' @param genes A character vector, input genes to be converted.
#' @param fromType The input ID type, one of "entrez", "symbol"(default), "hgnc",
#' "ensembl", "fullname" and "uniprotswissprot";
#' you can also input other valid attribute names for biomaRt.
#' Look at the code in examples to check valid attributes.
#' @param toType The output ID type, similar to `fromType`.
#' @param organism "hsa"(default), "mmu", "bta", "cfa", "ptr", "rno", and "ssc" are optional.
#' @param fromOrg "hsa", "mmu", "bta", "cfa", "ptr", "rno", and "ssc" are optional
#' (Only used when transform gene ids between organisms).
#' @param toOrg "hsa"(default), "mmu", "bta", "cfa", "ptr", "rno", and "ssc" are optional
#' (Only used when transform gene ids between organisms).
#' @param ensemblHost String, specifying ensembl host, you can use `listEnsemblArchives()`
#' to show all available Ensembl archives hosts.
#' @param update Boolean, specifying whether update built-in gene annotation (needs network and takes time).
#'
#' @return A character vector, named by unique input gene ids.
#'
#' @author Wubing Zhang
#'
#' @examples
#' data(mle.gene_summary)
#' TransGeneID(mle.gene_summary$Gene[1:10], organism="hsa")
#' TransGeneID(mle.gene_summary$Gene[1:10], toType="Symbol", fromOrg = "hsa", toOrg = "mmu")
#' @import biomaRt
#' @export

TransGeneID <- function(genes, fromType="Symbol", toType="Entrez",
                        organism = "hsa", fromOrg = organism, toOrg = organism,
                        ensemblHost = "www.ensembl.org", update = FALSE){

  #### Verify  parameters ####
  genes = as.character(genes)
  fromType = tolower(fromType)
  toType = tolower(toType)
  if(length(genes)<1) return(c())
  keggcode = rep(c("hsa", "mmu", "rno", "bta", "cfa", "ptr", "ssc"), 2)
  names(keggcode) = c(tolower(c("Human", "Mouse", "Rat", "Bovine", "Canine", "Chimp", "Pig")),
                      c("hsa", "mmu", "rno", "bta", "cfa", "ptr", "ssc"))
  if(!tolower(organism)%in%names(keggcode)) stop("Organism error ...")
  if(!tolower(fromOrg)%in%names(keggcode)) stop("fromOrg error ...")
  if(!tolower(toOrg)%in%names(keggcode)) stop("toOrg error ...")

  organism = keggcode[tolower(organism)]
  fromOrg = keggcode[tolower(fromOrg)]
  toOrg = keggcode[tolower(toOrg)]

  #### Read annotation file ####
  datasets = paste0(c("hsapiens", "mmusculus", "btaurus", "cfamiliaris",
                      "ptroglodytes", "rnorvegicus", "sscrofa"), "_gene_ensembl")
  if(fromOrg==toOrg){#### GeneID Transformation within organisms ####
    ## ID mapping.
    if(all(c(fromType, toType) %in% c("entrez", "symbol", "hgnc", "ensembl"))){
      ann <- getGeneAnn(organism, update)[, c(fromType, toType)]
    }else{
      requireNamespace("biomaRt")
      ds = datasets[grepl(organism, datasets)]
      ensembl <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = ds, host = ensemblHost)
      ## decide the attributes automatically
      attrs = listAttributes(ensembl)$name
      if(sum(attrs==fromType)==0){
        idx1 = grepl(tolower(fromType), attrs)
        idx = idx1
        if(sum(idx1)>2) idx = idx1&grepl("_id", attrs)
        fromType = ifelse(sum(idx)>0, attrs[idx][1], attrs[idx1][1])
        if(fromType=="hgnc_symbol" & fromOrg=="mmu") fromType = "mgi_symbol"
      }
      if(sum(attrs==toType)==0){
        idx1 = grepl(tolower(toType), attrs)
        idx = idx1
        if(sum(idx1)>2) idx = idx1&grepl("_id", attrs)
        toType = ifelse(sum(idx)>0, attrs[idx][1], attrs[idx1])
        if(toType=="hgnc_symbol" & toOrg=="mmu") toType = "mgi_symbol"
      }
      ## retrieve the data
      ann = getBM(attributes=c(fromType, toType), mart = ensembl,
                  filters = fromType, values = genes)
    }
    ## merge the annotation
    colnames(ann) = c(fromType, toType)
    ## Retain unique conversion
    idx = ann[, toType]=="" | is.na(ann[, toType])
    ann = ann[!idx, ]
    idx = duplicated(ann[, fromType])
    convert = ann[!idx, toType]
    names(convert) = ann[!idx, fromType]
    gene_after = as.character(convert[genes])
    names(gene_after) = genes
  }else{#### GeneID Transformation between organisms ####
    if(all(c(fromType, toType) %in% c("symbol", "entrez"))){
      ## read built-in annotation
      ann = getOrtAnn(fromOrg, toOrg, update)
      ann = ann[, c(paste0(fromOrg, "_", fromType), paste0(toOrg, "_", toType))]
      colnames(ann) = c(fromOrg, toOrg)
    }else{
      ## Ortholog ID mapping.
      from = useMart("ensembl", dataset = datasets[grepl(fromOrg, datasets)])
      to = useMart("ensembl", dataset = datasets[grepl(toOrg, datasets)])
      ## decide the attributes automatically
      attrs_1 = listAttributes(from)$name
      attrs_2 = listAttributes(to)$name
      if(sum(attrs_1==fromType)==0){
        idx1 = grepl(tolower(fromType), attrs_1)
        idx = idx1
        if(sum(idx1)>2) idx = idx1&grepl("_id", attrs_1)
        fromType = ifelse(sum(idx)>0, attrs_1[idx][1], attrs_1[idx1][1])
        if(fromType=="hgnc_symbol" & fromOrg=="mmu") fromType = "mgi_symbol"
      }
      if(sum(attrs_2==toType)==0){
        idx1 = grepl(tolower(toType), attrs_2)
        idx = idx1
        if(sum(idx1)>2) idx = idx1&grepl("_id", attrs_2)
        toType = ifelse(sum(idx)>0, attrs_2[idx][1], attrs_2[idx1])
        if(toType=="hgnc_symbol" & toOrg=="mmu") toType = "mgi_symbol"
      }
      ## retrieve the data
      ann = getLDS(attributes = fromType, mart = from,
                   filters = fromType, values = genes,
                   attributesL = toType, martL = to)
      colnames(ann) = c(fromOrg, toOrg)
    }

    ## Retain unique conversion
    idx = ann[, toOrg]=="" | is.na(ann[, toOrg])
    ann = ann[!idx, ]
    idx = duplicated(ann[, fromOrg])
    convert = ann[!idx, toOrg]
    names(convert) = ann[!idx, fromOrg]
    gene_after = as.character(convert[genes])
    names(gene_after) = genes
  }
  return(gene_after)
}

#' Retrieve gene annotations from the NCBI, HNSC, and Uniprot databases.
#'
#' @docType methods
#' @name getGeneAnn
#' @rdname getGeneAnn
#'
#' @param org Character, hsa (default), bta, cfa, mmu, ptr, rno, ssc are optional.
#' @param update Boolean, indicating whether download current annotation.
#' @return A data frame.
#'
#' @author Wubing Zhang
#'
#' @examples
#' ann = getGeneAnn("hsa")
#' head(ann)
#'
#' @export
#'
getGeneAnn <- function(org = "hsa", update = FALSE){
  #### Read rds file directly ####
  rdsann = file.path(system.file("extdata", package = "MAGeCKFlute"),
                     paste0("NCBI_HGNC_GeneID_Annotation_", org, ".rds"))
  if(file.exists(rdsann) & !update) return(readRDS(rdsann))

  #### NCBI gene annotation ####
  gzfile = c("Homo_sapiens.gene_info.gz", "Bos_taurus.gene_info.gz",
             "Canis_familiaris.gene_info.gz", "Mus_musculus.gene_info.gz",
             "Pan_troglodytes.gene_info.gz", "Rattus_norvegicus.gene_info.gz",
             "Sus_scrofa.gene_info.gz")
  names(gzfile) = c("hsa", "bta", "cfa", "mmu", "ptr", "rno", "ssc")
  locfname <- file.path(system.file("extdata", package = "MAGeCKFlute"), gzfile[org])
  if((!file.exists(locfname)) | update){
    ## Download gene information from NCBI ftp server
    refname <- paste0("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/", gzfile[org])
    download.file(refname, locfname, quiet = TRUE)
  }

  ## Reorder the mapping file
  ncbi_ann = read.csv(gzfile(locfname), sep = "\t", header = TRUE,
                  quote = "", stringsAsFactors = FALSE, comment.char = "")
  ncbi_ann = ncbi_ann[, c("GeneID", "Symbol", "Synonyms", "dbXrefs", "type_of_gene", "description")]
  colnames(ncbi_ann)[c(1,2,6)] = c("entrez", "symbol", "fullname")
  ncbi_ann$hgnc = gsub("\\|.*", "", gsub(".*HGNC:", "", ncbi_ann$dbXrefs))
  ncbi_ann$ensembl = gsub("\\|.*", "", gsub(".*Ensembl:", "", ncbi_ann$dbXrefs))
  ncbi_ann$hgnc[!grepl("HGNC", ncbi_ann$dbXrefs)] = ""
  ncbi_ann$ensembl[!grepl("Ensembl", ncbi_ann$dbXrefs)] = ""

  synonyms_row = matrix(unlist(apply(ncbi_ann, 1, function(x){
    tmp = unlist(strsplit(x[3], "[|]"))
    if(tmp[1]!="" & tmp[1]!="-") return(as.vector(rbind(x[1], tmp, x[7], x[8], x[6])))
    return(NULL)
  })) , ncol=5, byrow = TRUE)
  colnames(synonyms_row) = c("entrez", "symbol", "hgnc", "ensembl", "fullname")
  ncbi_ann = rbind(ncbi_ann[,c(1,2,7,8,6)], synonyms_row)
  ncbi_ann = ncbi_ann[,-5]

  #### HGNC gene annotation ####
  if(org=="hsa"){
    locfname2 = file.path(system.file("extdata", package = "MAGeCKFlute"), "HGNC_GeneID_annotation.txt.gz")
    if((!file.exists(locfname2)) | update){
      ## Download gene information from HGNC
      refname <- "ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt"
      download.file(refname, locfname2, quiet = TRUE)
    }
    ## Reorder the mapping file
    hgnc_ann = read.csv(gzfile(locfname2), sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE, comment.char = "")
    hgnc_ann = hgnc_ann[, c("entrez_id", "symbol", "hgnc_id", "ensembl_gene_id", "name",
                            "alias_symbol", "prev_symbol", "refseq_accession")]
    hgnc_ann$alias_symbol = paste0(hgnc_ann$alias_symbol, '|', hgnc_ann$prev_symbol)

    synonyms_row = matrix(unlist(apply(hgnc_ann, 1, function(x){
      tmp = unlist(strsplit(x[6], "\\|"))
      if(length(tmp)>0) return(as.vector(rbind(x[1], tmp, x[3], x[4])))
      return(NULL)
    })) , ncol=4, byrow = TRUE)
    colnames(synonyms_row) = c("entrez", "symbol", "hgnc", "ensembl")
    synonyms_row = synonyms_row[synonyms_row[,2]!="", ]
    names(hgnc_ann)[1:4] = c("entrez", "symbol", "hgnc", "ensembl")
    hgnc_ann$hgnc = gsub("HGNC:", "", hgnc_ann$hgnc)
    hgnc_ann = rbind(hgnc_ann[,1:4], synonyms_row)
  }

  #### Ensembl gene annotation ####
  datasets = paste0(c("hsapiens", "mmusculus", "btaurus", "cfamiliaris",
                      "ptroglodytes", "rnorvegicus", "sscrofa"), "_gene_ensembl")
  ds = datasets[grepl(org, datasets)]
  ensembl <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = ds)
  symbol <- ifelse(org=="mmu", "mgi_symbol", "hgnc_symbol")
  ensembl_ann = getBM(attributes=c("entrezgene_id", symbol, "hgnc_id", "ensembl_gene_id"), mart = ensembl)
  colnames(ensembl_ann) = c("entrez", "symbol", "hgnc", "ensembl")
  ensembl_ann$hgnc = gsub("HGNC:", "", ensembl_ann$hgnc)

  #### Merge HGNC and NCBI annotation ####
  if(org=="hsa"){
    data = rbind.data.frame(ncbi_ann, hgnc_ann, ensembl_ann)
  }else data = rbind.data.frame(ncbi_ann, ensembl_ann)
  data$entrez = as.character(as.integer(data$entrez))
  data$hgnc = gsub(" *", "", data$hgnc)
  idx = duplicated(paste(data$entrez, data$symbol, sep = "_"))
  data = data[!idx, ]
  rownames(data) = NULL
  saveRDS(data, rdsann)
  return(data)
}


#' Retreive reference orthologs annotation.
#'
#' @docType methods
#' @name getOrtAnn
#' @rdname getOrtAnn
#'
#' @param fromOrg Character, hsa (default), bta, cfa, mmu, ptr, rno, ssc are optional.
#' @param toOrg Character, hsa (default), bta, cfa, mmu, ptr, rno, ssc are optional.
#' @param update Boolean, indicating whether download recent annotation from NCBI.
#' @return A data frame.
#'
#' @author Wubing Zhang
#'
#' @examples
#' ann = getOrtAnn("hsa", "mmu")
#' head(ann)
#'
#' @export
#'
getOrtAnn <- function(fromOrg = "mmu", toOrg = "hsa", update = FALSE){
  #### Read rds file directly ####
  rdsann = file.path(system.file("extdata", package = "MAGeCKFlute"),
                     paste0("HOM_GeneID_Annotation_", fromOrg, "_", toOrg, ".rds"))
  if(file.exists(rdsann) & !update) return(readRDS(rdsann))

  keggcode = c("hsa", "mmu", "rno", "bta", "cfa", "ptr", "ssc")
  names(keggcode) = c("human", "mouse", "rat", "bovine", "canine", "chimp", "pig")
  #### Download data from MGI ####
  locfname <- file.path(system.file("extdata", package = "MAGeCKFlute"),
                        "HOM_MouseHumanSequence.rpt.gz")
  if((!file.exists(locfname)) | update){
    ## Download gene information
    refname <- "http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"
    download.file(refname, locfname, quiet = TRUE)
  }
  ## Reorder the mapping file
  read.table(locfname, sep = "\t", header = TRUE, stringsAsFactors = FALSE) -> mgi_ann
  mgi_ann = mgi_ann[, c(1,2,4,5)]
  colnames(mgi_ann) = c("homoloid", "org", "symbol", "entrez")
  mgi_ann$org = gsub(", laboratory", "", mgi_ann$org)
  mgi_ann$org = keggcode[mgi_ann$org]

  #### Download data from NCBI ####
  locfname <- file.path(system.file("extdata", package = "MAGeCKFlute"),
                        "homologene.data.gz")
  if((!file.exists(locfname)) | update){
    ## Download gene information
    refname <- "ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data"
    download.file(refname, locfname, quiet = TRUE)
  }
  ## Reorder the mapping file
  read.table(locfname, sep = "\t", stringsAsFactors = FALSE, quote = "") -> ncbi_ann
  ncbi_ann = ncbi_ann[, c(1,2,4,3)]
  colnames(ncbi_ann) = c("homoloid", "org", "symbol", "entrez")
  names(keggcode) = c(9606, 10090, 10116, 9913, 9615, 9598, 9823)
  ncbi_ann$org = keggcode[as.character(ncbi_ann$org)]
  ncbi_ann = ncbi_ann[!is.na(ncbi_ann$org), ]

  ## Merge and arrange the mapping file
  ann = rbind.data.frame(mgi_ann, ncbi_ann)
  genes = unique(ann$entrez)
  idx1 = ann$entrez %in% genes
  idx2 = ann$org == toOrg
  idx3 = ann$homoloid %in% ann$homoloid[idx1]
  tmp1 = ann[idx1, c("homoloid", "symbol", "entrez")]
  tmp2 = ann[(idx2&idx3), c("homoloid", "symbol", "entrez")]
  colnames(tmp1)[2:3] = paste0(fromOrg, c("_symbol", "_entrez"))
  colnames(tmp2)[2:3] = paste0(toOrg, c("_symbol", "_entrez"))
  ann = merge(tmp1, tmp2, by = "homoloid")[,-1]

  #### Retrieve annotation from Ensembl ####
  datasets = paste0(c("hsapiens", "mmusculus", "btaurus", "cfamiliaris",
                      "ptroglodytes", "rnorvegicus", "sscrofa"), "_gene_ensembl")
  ## Ortholog ID mapping.
  from = useMart("ensembl", dataset = datasets[grepl(fromOrg, datasets)])
  to = useMart("ensembl", dataset = datasets[grepl(toOrg, datasets)])
  ## decide the attributes automatically
  from_symbol <- ifelse(fromOrg=="mmu", "mgi_symbol", "hgnc_symbol")
  to_symbol <- ifelse(toOrg=="mmu", "mgi_symbol", "hgnc_symbol")
  ## retrieve the data
  ensembl_ann = getLDS(attributes = c(from_symbol, "entrezgene_id"), mart = from,
                       attributesL = c(to_symbol, "entrezgene_id"), martL = to)
  colnames(ensembl_ann) = c(paste0(fromOrg, c("_symbol", "_entrez")),
                            paste0(toOrg, c("_symbol", "_entrez")))

  ## Merge all the annotations
  ann = rbind.data.frame(ann, ensembl_ann)
  idx = duplicated(paste(ann[,1], ann[,2], ann[,3], ann[,4], sep = "_"))
  ann = ann[!idx, ]
  saveRDS(ann, rdsann)
  return(ann)
}




#' Get the kegg code of specific mammalia organism.
#'
#' @docType methods
#' @name getOrg
#' @rdname getOrg
#'
#' @param organism Character, KEGG species code, or the common species name.
#' For all potential values check: data(bods); bods. Default org="hsa",
#' and can also be "human" (case insensitive).
#' @return A list containing three elements:
#' \item{org}{species}
#' \code{pkg}{annotation package name}
#'
#' @author Wubing Zhang
#'
#' @examples
#' ann = getOrg("human")
#' print(ann$pkg)
#'
#' @export

getOrg <- function(organism){
  requireNamespace("clusterProfiler", quietly=TRUE) || stop("need clusterProfiler package")
  bods <- data.frame(package = paste0("org.", c("Hs", "Mm", "Rn", "Bt", "Cf", "Pt", "Ss"), ".eg.db"),
                     species = c("Human", "Mouse", "Rat", "Bovine", "Canine", "Chimp", "Pig"),
                     "kegg code" = c("hsa", "mmu", "rno", "bta", "cfa", "ptr", "ssc"),
                     check.names = FALSE, stringsAsFactors = FALSE)
  res=list()
  ##======
  # Get the mapping from organism to package
  ridx = c(which(tolower(bods[,2])==tolower(organism)),
           which(tolower(bods[,3])==tolower(organism)))
  stopifnot(length(ridx)==1)
  res$org = bods[ridx,3]
  res$pkg = bods[ridx,1]
  return(res)
}
