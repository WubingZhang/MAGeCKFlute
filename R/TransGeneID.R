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
#' @param unique Boolean, specifying whether do one-to-one mapping.
#' @param update Boolean, specifying whether update built-in gene annotation (needs network and takes time).
#'
#' @return A character vector, named by unique input gene ids.
#'
#' @author Wubing Zhang
#'
#' @examples
#' TransGeneID("HLA-A", organism="hsa")
#' TransGeneID("HLA-A", toType = "uniprot", organism="hsa")
#' TransGeneID("H2-K1", toType="Symbol", fromOrg = "mmu", toOrg = "hsa")
#'
#' @export

TransGeneID <- function(genes, fromType="Symbol", toType="Entrez",
                        organism = "hsa", fromOrg = organism, toOrg = organism,
                        ensemblHost = "www.ensembl.org",
                        unique = TRUE, update = FALSE){

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
    if(all(c(fromType, toType) %in% c("entrez", "symbol", "ensembl"))){
      ann <- getGeneAnn(organism, update=update)$Gene
      if("symbol" %in% c(fromType, toType) & any(!genes%in%ann$symbol)){
        ann = rbind(ann, ann)
        idx = (nrow(ann)/2+1):nrow(ann)
        ann$symbol[idx] = ann$synonyms[idx]
      }
      ann = ann[, c(fromType, toType)]
    }else if(all(c(fromType, toType) %in% c("uniprot", "ensembl", "refseq", "symbol"))){
      ann <- getGeneAnn(organism, update=update)$Protein
      ann = ann[, c(fromType, toType)]
    }else{
      if (!requireNamespace("biomaRt", quietly = TRUE)) {
        stop("Package \"biomaRt\" is required. Please install it.", call. = FALSE)
      }
      ds = datasets[grepl(organism, datasets)]
      ensembl <- biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL',
                                  dataset = ds, host = ensemblHost)
      ## decide the attributes automatically
      attrs = biomaRt::listAttributes(ensembl)$name
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
      ann = biomaRt::getBM(attributes=c(fromType, toType), mart = ensembl,
                           filters = fromType, values = genes)
    }

    ## Retain unique conversion
    idx = ann[, toType]=="" | is.na(ann[, toType])
    ann = ann[!idx, ]
    idx = ann[, fromType]=="" | is.na(ann[, fromType])
    ann = ann[!idx, ]
    ##
    tmp = ann
    tmp[, fromType] = gsub("\\..*|-.*", "", tmp[, fromType])
    ann = rbind.data.frame(ann, tmp)
    ann = ann[ann[,fromType]%in%genes, ]
    # idx = duplicated(ann[, fromType])
    # convert = ann[!idx, toType]
    # names(convert) = ann[!idx, fromType]
    # gene_after = as.character(convert[genes])
    # names(gene_after) = genes
    # genes = gsub("\\..*|-.*", "", genes)
    # gene_after[is.na(gene_after)] = as.character(convert[genes[is.na(gene_after)]])
  }else{#### GeneID Transformation between organisms ####
    if(all(c(fromType, toType) %in% c("symbol", "entrez"))){
      ## read built-in annotation
      ann = getOrtAnn(fromOrg, toOrg, update)
      ann = ann[, c(paste0(fromOrg, "_", fromType), paste0(toOrg, "_", toType))]
      colnames(ann) = c(fromOrg, toOrg)
    }else{
      if (!requireNamespace("biomaRt", quietly = TRUE)) {
        stop("Package \"biomaRt\" is required. Please install it.",
             call. = FALSE)
      }
      ## Ortholog ID mapping.
      from = biomaRt::useMart("ensembl", dataset = datasets[grepl(fromOrg, datasets)])
      to = biomaRt::useMart("ensembl", dataset = datasets[grepl(toOrg, datasets)])
      ## decide the attributes automatically
      attrs_1 = biomaRt::listAttributes(from)$name
      attrs_2 = biomaRt::listAttributes(to)$name
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
      ann = biomaRt::getLDS(attributes = fromType, mart = from,
                            filters = fromType, values = genes,
                            attributesL = toType, martL = to)
      colnames(ann) = c(fromOrg, toOrg)
    }

    ## Retain unique conversion
    idx = ann[, toOrg]=="" | is.na(ann[, toOrg])
    ann = ann[!idx, ]
    ann = ann[ann[,fromOrg]%in%genes, ]
    # idx = duplicated(ann[, fromOrg])
    # convert = ann[!idx, toOrg]
    # names(convert) = ann[!idx, fromOrg]
    # gene_after = as.character(convert[genes])
    # names(gene_after) = genes
  }
  ann = ann[!duplicated(paste0(ann[,1], ann[,2])), ]
  if(unique){
    ann = ann[!duplicated(ann[,1]), ]
    rownames(ann) = ann[,1]
    tmp = ann[genes,2]
    names(tmp) = genes
    ann = tmp
  }
  return(ann)
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
#' \dontrun{
#'   ann = getGeneAnn("hsa")
#'   head(ann)
#' }
#'
#' @export
#'
getGeneAnn <- function(org = "hsa", update = FALSE){
  options(stringsAsFactors = FALSE)
  #### Read rds file directly ####
  rdsann = file.path(system.file("extdata", package = "MAGeCKFlute"),
                     paste0("GeneID_Annotation_", org, ".rds"))
  if(file.exists(rdsann) & !update)
    return(list(Gene = readRDS(rdsann), Protein = readRDS(gsub("Gene", "Protein", rdsann))))

  #### NCBI gene annotation ####
  gzfile = paste0(c("Homo_sapiens", "Bos_taurus", "Canis_familiaris", "Mus_musculus",
                    "Pan_troglodytes", "Rattus_norvegicus", "Sus_scrofa"), ".gene_info.gz")
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

  ncbi_ann = matrix(unlist(apply(ncbi_ann, 1, function(x){
    tmp = unlist(strsplit(x[3], "[|]"))
    return(as.vector(rbind(x[1], x[2], tmp, x[7], x[8], x[6])))
  })) , ncol=6, byrow = TRUE)
  colnames(ncbi_ann) = c("entrez", "symbol", "synonyms", "hgnc", "ensembl", "fullname")
  ncbi_ann[,1] = gsub(" ", "", ncbi_ann[,1])

  #### HGNC gene annotation ####
  # if(org=="hsa"){
  #   locfname2 = file.path(system.file("extdata", package = "MAGeCKFlute"), "HGNC_GeneID_annotation.txt.gz")
  #   if((!file.exists(locfname2)) | update){
  #     ## Download gene information from HGNC
  #     refname <- "ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt"
  #     download.file(refname, locfname2, quiet = TRUE)
  #   }
  #   ## Reorder the mapping file
  #   hgnc_ann = read.csv(gzfile(locfname2), sep = "\t", header = TRUE,
  #                       stringsAsFactors = FALSE, comment.char = "")
  #   hgnc_ann = hgnc_ann[, c("entrez_id", "ensembl_gene_id", "symbol", "hgnc_id", "name",
  #                           "alias_symbol", "prev_symbol", "refseq_accession")]
  #   hgnc_ann$alias_symbol = paste0(hgnc_ann$alias_symbol, '|', hgnc_ann$prev_symbol)
  #
  #   synonyms_row = matrix(unlist(apply(hgnc_ann, 1, function(x){
  #     tmp = unlist(strsplit(x[6], "\\|"))
  #     if(length(tmp)>0) return(as.vector(rbind(x[1], tmp, x[2], x[8])))
  #     return(NULL)
  #   })) , ncol=4, byrow = TRUE)
  #   colnames(synonyms_row) = c("entrez", "symbol", "ensembl", "refseq")
  #   synonyms_row = synonyms_row[synonyms_row[,2]!="", ]
  #   names(hgnc_ann)[1:3] = c("entrez", "ensembl", "symbol")
  #   names(hgnc_ann)[8] = "refseq"
  #   hgnc_ann = rbind(hgnc_ann[,c(1:3,8)], synonyms_row[,c(1,3,2,4)])
  #   ncbi_ann = merge(ncbi_ann, hgnc_ann, by = names(hgnc_ann)[1:3], all = TRUE)
  # }

  #### Ensembl gene annotation ####
  tmpfile = file.path(system.file("extdata", package = "MAGeCKFlute"), "filelist")
  download.file("ftp://ftp.ensembl.org/pub/", tmpfile, quiet = TRUE)
  tmp = read.table(tmpfile, fill = TRUE, quote = "", stringsAsFactors = FALSE)
  tmp = gsub("release-", "", tmp[grepl("release", tmp[, ncol(tmp)]), ncol(tmp)])
  version = tmp[length(tmp)]
  gzfile = c("Homo_sapiens.GRCh38.", "Bos_taurus.ARS-UCD1.2.", "Canis_familiaris.CanFam3.1.",
             "Mus_musculus.GRCm38.", "Pan_troglodytes.Pan_tro_3.0.",
             "Rattus_norvegicus.Rnor_6.0.", "Sus_scrofa.Sscrofa11.1.")
  names(gzfile) = c("hsa", "bta", "cfa", "mmu", "ptr", "rno", "ssc")
  entrezfile <- paste0("ftp://ftp.ensembl.org/pub/release-", version, "/tsv/",
                       tolower(gsub("\\..*", "", gzfile[org])), "/", gzfile[org],
                       version, ".entrez.tsv.gz")
  # uniprotfile <- paste0("ftp://ftp.ensembl.org/pub/release-", version, "/tsv/",
  #                       tolower(gsub("\\..*", "", gzfile[org])), "/", gzfile[org],
  #                       version, ".uniprot.tsv.gz")
  refseqfile <- paste0("ftp://ftp.ensembl.org/pub/release-", version, "/tsv/",
                       tolower(gsub("\\..*", "", gzfile[org])), "/", gzfile[org],
                       version, ".refseq.tsv.gz")
  download.file(entrezfile, tmpfile, quiet = TRUE)
  ensg_entrez = read.table(tmpfile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  ensg_entrez = ensg_entrez[, c(1,4)]
  ensg_entrez = ensg_entrez[!duplicated(paste0(ensg_entrez[,1], ensg_entrez[,2])), ]
  # download.file(uniprotfile, tmpfile, quiet = TRUE)
  # ensg_uniprot = read.table(tmpfile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  # ensg_uniprot = ensg_uniprot[, c(1,4)]
  # ensg_uniprot = ensg_uniprot[!duplicated(paste0(ensg_uniprot[,1], ensg_uniprot[,2])), ]
  download.file(refseqfile, tmpfile, quiet = TRUE)
  ensg_refseq = read.table(tmpfile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  ensg_refseq = ensg_refseq[, c(1,4)]
  ensg_refseq = ensg_refseq[!duplicated(paste0(ensg_refseq[,1], ensg_refseq[,2])), ]
  ensembl_ann = merge(ensg_entrez, ensg_refseq, by = names(ensg_entrez)[1], all = TRUE)
  colnames(ensembl_ann) = c("ensembl", "entrez", "refseq")
  ## Remove redundant refseq ids.
  idx = duplicated(paste(ensembl_ann$ensembl, ensembl_ann$entrez))
  ensembl_ann = ensembl_ann[!idx, ]
  idx = is.na(ensembl_ann$entrez) & is.na(ensembl_ann$refseq)
  ensembl_ann = ensembl_ann[!idx, ]
  ensembl_ann = ensembl_ann[grepl("ENSG", ensembl_ann$ensembl), ]
  # datasets = paste0(c("hsapiens", "mmusculus", "btaurus", "cfamiliaris",
  #                     "ptroglodytes", "rnorvegicus", "sscrofa"), "_gene_ensembl")
  # ds = datasets[grepl(org, datasets)]
  # ensembl <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = ds)
  # symbol <- ifelse(org=="mmu", "mgi_symbol", "hgnc_symbol")
  # ensembl_ann = getBM(attributes=c("entrezgene_id", symbol, "hgnc_id", "ensembl_gene_id"), mart = ensembl)
  # colnames(ensembl_ann) = c("entrez", "symbol", "hgnc", "ensembl")
  # ensembl_ann$hgnc = gsub("HGNC:", "", ensembl_ann$hgnc)

  #### Merge all annotations ####
  geneann = merge(ncbi_ann[, c("ensembl","entrez","symbol","synonyms")],
                  ensembl_ann[,c("ensembl","entrez")], by = c("ensembl","entrez"), all = TRUE)
  geneann$entrez = gsub(" ", "", geneann$entrez)
  # ids = gsub("-.*$", "", uniprot_ann$Canonical[!grepl("-1", uniprot_ann$Canonical)])
  # data$uniprot[data$uniprot%in%ids] = uniprot_ann[data$uniprot[data$uniprot%in%ids], "Canonical"]
  # idx = !(grepl("-", data$uniprot)|is.na(data$uniprot))
  # data$uniprot[idx] = paste0(data$uniprot[idx], "-1")
  saveRDS(geneann, rdsann)

  #### Uniprot gene annotation ####
  proteome_code = c("up000005640", "UP000009136", "UP000002254", "up000000589", "UP000002277", "UP000002494")
  names(proteome_code) = c("hsa", "bta", "cfa", "mmu", "ptr", "rno")
  locfname <- file.path(system.file("extdata", package = "MAGeCKFlute"),
                        paste0("uniprot_proteome_", proteome_code[org], ".tab"))
  uniprot_link <- paste0("https://www.uniprot.org/uniprot/?query=proteome:", proteome_code[org],
                         "&format=tab&force=true&columns=id,entry%20name,reviewed,protein%20names,genes,organism,database(Ensembl),comment(SUBCELLULAR%20LOCATION),database(RefSeq),comment(ALTERNATIVE%20PRODUCTS)&sort=score")

  if((!file.exists(locfname)) | update){
    ## Download protein annotation information from uniprot
    download.file(uniprot_link, locfname, quiet = TRUE)
  }
  ## Reorder the mapping file
  uniprot_ann = read.csv(gzfile(locfname), sep = "\t", header = TRUE,
                         quote = "", stringsAsFactors = FALSE, comment.char = "")
  suppressWarnings(try(file.remove(locfname), silent = TRUE))
  colnames(uniprot_ann) = c("Entry", "EntryName", "Status", "Name", "Gene",
                            "Organism", "Ensembl", "Subcellular", "RefSeq", "Isoforms")
  uniprot_ann$Canonical = gsub(".*IsoId=", "", gsub("; Sequence=Displayed.*", "", uniprot_ann$Isoforms))
  uniprot_ann$Canonical[!grepl("Sequence=Displayed", uniprot_ann$Isoforms)] =
    paste0(uniprot_ann$Entry[!grepl("Sequence=Displayed", uniprot_ann$Isoforms)], "-1")
  rownames(uniprot_ann) = uniprot_ann$Entry
  Symbols = unlist(apply(uniprot_ann, 1, function(x){
    Uniprot = x[1]
    Symbols = unlist(strsplit(gsub(";$", "", x[5]), " "))
    Symbols = Symbols[Symbols!=""|Symbols==""]
    if(length(Symbols)==0) return(NULL)
    if(length(Symbols)==1) Symbols = c(Symbols, "")
    rbind(rep(Uniprot, length(Symbols)-1), Symbols[1], Symbols[-1])
  }))
  Symbols = matrix(Symbols, ncol = 3, byrow = TRUE)
  colnames(Symbols) = c('Entry', "symbol", "synonyms")
  Symbols = as.data.frame(Symbols, stringsAsFactors = FALSE)

  ENSTs = unlist(apply(uniprot_ann, 1, function(x){
    Uniprot = x[1]
    ENSTs = unlist(strsplit(gsub(";$", "", x[7]), ";"))
    rbind(rep(Uniprot, length(ENSTs)), ENSTs)
  }))
  ENSTs = matrix(ENSTs, ncol = 2, byrow = TRUE)
  ENSTs = as.data.frame(ENSTs, stringsAsFactors = FALSE)
  colnames(ENSTs) = c('Entry', "ensembl")
  ENSTs$uniprot = gsub(".*\\[|\\]", "", ENSTs$ensembl)
  ENSTs$uniprot[!grepl("\\[", ENSTs$ensembl)] = ENSTs$Entry[!grepl("\\[", ENSTs$ensembl)]
  ENSTs$ensembl = gsub(" .*", "", ENSTs$ensembl)
  ENSTs$uniprot[!grepl("-", ENSTs$uniprot)] =
    uniprot_ann[ENSTs$uniprot[!grepl("-", ENSTs$uniprot)], "Canonical"]

  RefSeq = unlist(apply(uniprot_ann, 1, function(x){
    Uniprot = x[1]
    RefSeq = unlist(strsplit(gsub(";$", "", x[9]), ";"))
    rbind(rep(Uniprot, length(RefSeq)), RefSeq)
  }))
  RefSeq = matrix(RefSeq, ncol = 2, byrow = TRUE)
  RefSeq = as.data.frame(RefSeq, stringsAsFactors = FALSE)
  colnames(RefSeq) = c('Entry', "refseq")
  RefSeq$uniprot = gsub(".*\\[|\\]", "", RefSeq$refseq)
  RefSeq$uniprot[!grepl("\\[", RefSeq$refseq)] = RefSeq$Entry[!grepl("\\[", RefSeq$refseq)]
  RefSeq$refseq = gsub(" .*", "", RefSeq$refseq)
  RefSeq$uniprot[!grepl("-", RefSeq$uniprot)] =
    uniprot_ann[RefSeq$uniprot[!grepl("-", RefSeq$uniprot)], "Canonical"]
  ## Merge all annotations from each column
  uniprot_ann = merge(ENSTs[,-1], RefSeq[,-1], by = "uniprot", all = TRUE)
  uniprot_ann$Entry = gsub("-.*","",uniprot_ann$uniprot)
  # Symbols = Symbols[!duplicated(Symbols$Entry), ]
  uniprot_ann = merge(uniprot_ann, Symbols, by = "Entry", all = TRUE)
  saveRDS(uniprot_ann, gsub("Gene", "Protein", rdsann))

  return(list(Gene = geneann, Protein = uniprot_ann))
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
#' \dontrun{
#'   ann = getOrtAnn("mmu", "hsa")
#'   head(ann)
#' }
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
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Package \"biomaRt\" is required. Please install it.", call. = FALSE)
  }
  from = biomaRt::useMart("ensembl", dataset = datasets[grepl(fromOrg, datasets)])
  to = biomaRt::useMart("ensembl", dataset = datasets[grepl(toOrg, datasets)])
  ## decide the attributes automatically
  from_symbol <- ifelse(fromOrg=="mmu", "mgi_symbol", "hgnc_symbol")
  to_symbol <- ifelse(toOrg=="mmu", "mgi_symbol", "hgnc_symbol")
  ## retrieve the data
  ensembl_ann = biomaRt::getLDS(attributes = c(from_symbol, "entrezgene_id"), mart = from,
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
