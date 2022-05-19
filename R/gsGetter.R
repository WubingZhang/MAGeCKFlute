#' Extract pathway annotation from GMT file.
#'
#' @docType methods
#' @name gsGetter
#' @rdname gsGetter
#'
#' @param gmtpath The path to customized gmt file.
#' @param type Molecular signatures for testing, available datasets include
#' Pathway (KEGG, REACTOME, C2_CP:PID, C2_CP:BIOCARTA), GO (GOBP, GOCC, GOMF),
#' MSIGDB (C1, C2 (C2_CP (C2_CP:PID, C2_CP:BIOCARTA), C2_CGP),
#' C3 (C3_MIR, C3_TFT), C4 (C4_CGN, C4_CM), C5 (C5_BP, C5_CC, C5_MF), C6, C7, H)
#' and Complex (CORUM). Any combination of them are also accessible
#' (e.g. 'GOBP+GOMF+KEGG+REACTOME').
#' @param limit A two-length vector, specifying the minimal and
#' maximal size of gene sets to load.
#' @param organism 'hsa' or 'mmu'.
#' @param update Boolean, indicating whether update the gene sets from source database.
#'
#' @return A three-column data frame.
#'
#' @author Wubing Zhang
#'
#' @examples
#' gene2path = gsGetter(type = "REACTOME+KEGG")
#' head(gene2path)
#'
#' @export
#'
gsGetter <- function(gmtpath = NULL, type = "All", limit = c(0, Inf),
                     organism = 'hsa', update = FALSE){
  ## Normalize type
  type = toupper(unlist(strsplit(type, "\\+")))
  if("ALL" %in% type) type = c("PATHWAY", "GO", "COMPLEX", "MSIGDB")
  if("MSIGDB" %in% type) type = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "H", type)
  if("GO" %in% type) type = c("GOBP", "GOCC", "GOMF", type)
  if("C2" %in% type) type = c("KEGG", "REACTOME", "C2_CP:PID", "C2_CP:BIOCARTA", "C2_CGP", type)
  if("PATHWAY" %in% type) type = c("KEGG", "REACTOME", "C2_CP:PID", "C2_CP:BIOCARTA", type)
  if("COMPLEX" %in% type) type = c("CORUM", type)
  type = setdiff(type, c("ALL", "MSIGDB", "C2", "GO", "PATHWAY", "COMPLEX"))

  ## Update genesets
  if(update) retrieve_gs(organism=organism, type = type)

  ## read GMT files
  if(!is.null(gmtpath)){
    gene2path = ReadGMT(gmtpath, limit = limit)
  }else{
    gene2path = data.frame()
    if("KEGG" %in% type){
      gsfile = file.path(system.file("extdata", package = "MAGeCKFlute"),
                         paste0("kegg.all.entrez.", organism, ".rds"))
      if(!file.exists(gsfile)) retrieve_gs(type = "KEGG", organism=organism)
      tmp = readRDS(gsfile)
      colnames(tmp) = c("ENTREZID", "PathwayID", "PathwayName")
      gene2path = rbind(gene2path, tmp)
    }
    if("REACTOME" %in% type){
      gsfile = file.path(system.file("extdata", package = "MAGeCKFlute"),
                         paste0("reactome.all.entrez.", organism, ".rds"))
      if(!file.exists(gsfile)) retrieve_gs(type = "REACTOME", organism=organism)
      tmp = readRDS(gsfile)
      colnames(tmp) = c("ENTREZID", "PathwayID", "PathwayName")
      gene2path = rbind(gene2path, tmp)
    }
    if("CORUM" %in% type){
      gsfile = file.path(system.file("extdata", package = "MAGeCKFlute"),
                         paste0("corum.all.entrez.", organism, ".rds"))
      if(!file.exists(gsfile)) retrieve_gs(type = "CORUM", organism=organism)
      tmp = readRDS(gsfile)
      colnames(tmp) = c("ENTREZID", "PathwayID", "PathwayName")
      gene2path = rbind(gene2path, tmp)
    }
    if(any(grepl("^GO", type))){
      gsfile = file.path(system.file("extdata", package = "MAGeCKFlute"),
                         paste0("go.all.entrez.", organism, ".rds"))
      if(!file.exists(gsfile)) retrieve_gs(type = "GO", organism=organism)
      go = readRDS(gsfile)
      go = go[go$Category%in%gsub("GO", "", type), 1:3]
      colnames(go) = c("ENTREZID", "PathwayID", "PathwayName")
      gene2path = rbind(gene2path, go)
    }
    othertypes = setdiff(type, c("KEGG", "CORUM", "REACTOME", "GOBP", "GOMF", "GOCC"))
    if(length(othertypes)>0){
      if (!requireNamespace("msigdbr", quietly = TRUE)) {
        stop("Package \"msigdbr\" is required. Please install it.", call. = FALSE)
      }
      for(i in othertypes){
        category = gsub("_.*", "", i)
        subcat  = NULL
        if(grepl("_", i)) subcat = gsub(".*_", "", i)
        species = ifelse(organism=="mmu", "Mus musculus", "Homo sapiens")
        m = msigdbr::msigdbr(species = species, category = category,
                             subcategory = subcat)[, c("entrez_gene", "gs_name", "gs_name")]
        colnames(m) = c("ENTREZID", "PathwayID", "PathwayName")
        m$PathwayName = gsub("^[[:alnum:]]*_", "", m$PathwayName)
        m$PathwayName = gsub("_", " ", m$PathwayName)
        # m$PathwayName = paste0(substr(m$PathwayName,0,1),
        #                        tolower(substr(m$PathwayName,2,nchar(m$PathwayName))))
        gene2path = rbind(gene2path, m)
      }
    }
    gene2path = na.omit(gene2path)
    count_gene = table(gene2path$PathwayID)
    pathways = names(count_gene)[count_gene>limit[1] & count_gene<=limit[2]]
    gene2path = gene2path[gene2path$PathwayID%in%pathways, ]
  }
  names(gene2path) = c("Gene","PathwayID", "PathwayName")
  # gene2path$PathwayName = tolower(gsub("_", " ", gene2path$PathwayName))
  gene2path$Gene = as.character(gene2path$Gene)
  return(gene2path)
}
#' Update genesets from source database
#'
#' @docType methods
#' @name retrieve_gs
#' @rdname retrieve_gs
#'
#' @param type A vector of databases, such as KEGG, REACTOME, CORUM, GO.
#' @param organism 'hsa' or 'mmu'.
#'
#' @return save data to local library.
#'
#' @author Wubing Zhang
#'
#' @export
#'
retrieve_gs <- function(type = c("KEGG", "REACTOME", "CORUM", "GO"), organism = 'hsa'){
  options(stringsAsFactors = FALSE)

  if("KEGG" %in% type){ ## Process genesets from KEGG
    message(format(Sys.time(), " Downloading genesets from KEGG ..."))
    gene2path = read.table(paste0("http://rest.kegg.jp/link/pathway/", organism),
                           sep = "\t", stringsAsFactors = FALSE)
    names(gene2path) = c("EntrezID", "PathwayID")
    pathways = read.table(paste0("http://rest.kegg.jp/list/pathway/", organism),
                          sep = "\t", stringsAsFactors = FALSE)
    names(pathways) = c("PathwayID","PathwayName")
    gene2path$EntrezID=gsub(".*:", "", gene2path$EntrezID)
    gene2path$PathwayID=gsub(".*:","",gene2path$PathwayID)
    pathways$PathwayID=gsub(".*:","",pathways$PathwayID)
    pathways$PathwayName=gsub(" - .*", "", pathways$PathwayName)
    rownames(pathways) = pathways$PathwayID
    gene2path$PathwayName = pathways[gene2path$PathwayID, "PathwayName"]
    locfname = file.path(system.file("extdata", package = "MAGeCKFlute"),
                         paste0("kegg.all.entrez.", organism, ".rds"))
    gene2path$PathwayID = paste0("KEGG_", gene2path$PathwayID)
    saveRDS(gene2path, locfname)
  }
  if("CORUM" %in% type){ ## Process genesets from CORUM
    message(format(Sys.time(), " Downloading genesets from CORUM ..."))
    base_url = "https://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt.zip"
    locfname = file.path(system.file("extdata", package = "MAGeCKFlute"), "allComplexes.txt.zip")
    download.file(base_url, locfname, quiet = TRUE)
    corum <- read.table(unz(locfname, "allComplexes.txt"), sep = "\t",
                        header = TRUE, quote = "", stringsAsFactors = FALSE)
    file.remove(locfname)
    if(organism=="hsa"){
      corum = corum[corum$Organism=="Human", ]
    }else if(organism=="mmu"){
      corum = corum[corum$Organism=="Mouse", ]
    }
    genes = strsplit(corum$subunits.Gene.name., ";")
    nset = unlist(lapply(genes, length))
    gene2corum = data.frame(EntrezID = unlist(genes),
                            ComplexID = rep(corum$ComplexID, nset),
                            ComplexName = rep(corum$ComplexName, nset))
    gene2corum$ComplexID = paste0("CORUM_", gene2corum$ComplexID)
    gene2corum$EntrezID = TransGeneID(gene2corum$EntrezID, "Symbol",
                                      "Entrez", organism = organism)
    gene2corum = na.omit(gene2corum)
    locfname = file.path(system.file("extdata", package = "MAGeCKFlute"),
                         paste0("corum.all.entrez.", organism, ".rds"))
    saveRDS(gene2corum, locfname)
  }
  if("REACTOME" %in% type){ ## Process genesets from REACTOME
    message(format(Sys.time(), " Downloading genesets from REACTOME ..."))
    gene2path = read.table("https://reactome.org/download/current/NCBI2Reactome.txt",
                           sep = "\t", stringsAsFactors = FALSE, comment.char = "", quote = "")
    colnames(gene2path) = c("EntrezID", "PathwayID", "link", "PathwayName", "Evidence", "Organism")
    gene2path = gene2path[grepl(organism, gene2path$PathwayID, ignore.case = TRUE), ]
    gene2path = gene2path[, c(1,2,4)]
    gene2path$PathwayID = gsub(paste0("R-", toupper(organism), "-"), "REACTOME_", gene2path$PathwayID)
    locfname = file.path(system.file("extdata", package = "MAGeCKFlute"),
                         paste0("reactome.all.entrez.", organism, ".rds"))
    saveRDS(gene2path, locfname)
  }
  ## Process genesets from Gene ontology
  if(any(grepl("^GO", type))){
    message(format(Sys.time(), " Downloading genesets from Gene Ontology ..."))
    tmpfile = file.path(system.file("extdata", package = "MAGeCKFlute"), "gene2go.gz")
    download.file("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz", destfile = tmpfile, quiet = TRUE)
    go <- read.table(gzfile(tmpfile), sep = "\t", header = TRUE,
                     stringsAsFactors = FALSE, comment.char = "", quote = "")
    colnames(go) <- c("tax_id","EntrezID","GO_ID","Evidence","Qualifier","GO_term", "PubMed","Category")
    taxid = c("hsa"=9606, "mmu"=10090)
    go = go[go$tax_id==taxid[organism], ]
    go = unique(go[, c(2,3,6,8)])
    go$Category[go$Category=="Process"] = "BP"
    go$Category[go$Category=="Component"] = "CC"
    go$Category[go$Category=="Function"] = "MF"
    go$EntrezID = as.character(go$EntrezID)
    locfname = file.path(system.file("extdata", package = "MAGeCKFlute"),
                         paste0("go.all.entrez.", organism, ".rds"))
    saveRDS(go, locfname)
  }
}
