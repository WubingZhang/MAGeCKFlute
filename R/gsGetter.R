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
#' gene2path = gsGetter(type = "REACTOME+CORUM")
#' head(gene2path)
#'
#' @importFrom msigdbr msigdbr
#' @export
#'
gsGetter <- function(gmtpath = NULL, type = "All", limit = c(0, Inf),
                     organism = 'hsa', update = FALSE){
  ## Update genesets
  if(update) retrieve_gs(organism=organism)
  ## Normalize type
  type = toupper(unlist(strsplit(type, "\\+")))
  if("ALL" %in% type) type = c("PATHWAY", "GO", "COMPLEX", "MSIGDB")
  if("MSIGDB" %in% type) type = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "H", type)
  if("GO" %in% type) type = c("C5", type)
  if("C2" %in% type) type = c("KEGG", "REACTOME", "C2_CP:PID", "C2_CP:BIOCARTA", "C2_CGP", type)
  if("PATHWAY" %in% type) type = c("KEGG", "REACTOME", "C2_CP:PID", "C2_CP:BIOCARTA", type)
  if("COMPLEX" %in% type) type = c("CORUM", type)
  type = setdiff(type, c("ALL", "MSIGDB", "C2", "GO", "PATHWAY", "COMPLEX"))
  type = gsub("GO", "C5_", type)
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
    if(length(setdiff(type, c("KEGG", "CORUM", "REACTOME")))>0){
      for(i in setdiff(type, c("KEGG", "CORUM", "REACTOME"))){
        category = gsub("_.*", "", i)
        subcat = gsub(".*_", "", i)
        if(subcat=="") subcat = NULL
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
#' @param type A vector of databases, such as KEGG, REACTOME, CORUM.
#' @param organism 'hsa' or 'mmu'.
#'
#' @return save data to local library.
#'
#' @author Wubing Zhang
#'
#' @export
#'
retrieve_gs <- function(type = c("KEGG", "REACTOME", "CORUM"), organism = 'hsa'){
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
    base_url = "http://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt.zip"
    locfname = file.path(system.file("extdata", package = "MAGeCKFlute"), "allComplexes.txt.zip")
    download.file(base_url, locfname, quiet = TRUE)
    corum <- read.table(unz(locfname, "allComplexes.txt"), sep = "\t",
                        header = TRUE, quote = "", stringsAsFactors = FALSE)
    file.remove(locfname)
    if(organism=="hsa")
      corum = corum[corum$Organism=="Human", ]
    else if(organism=="mmu")
      corum = corum[corum$Organism=="Mouse", ]
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
#   ## Molecular signature database
#   message(format(Sys.time(), " Downloading genesets from MsigDB ..."))
#   msigfile = file.path(system.file("extdata", package = "MAGeCKFlute"),
#                          paste0(organism, "_msig_entrez.gmt.gz"))
#   gene2path = ReadGMT(msigfile, limit = c(0, Inf))
#   saveRDS(gene2path, file.path(system.file("extdata", package = "MAGeCKFlute"),
#                                paste0("msigdb.all.entrez.", organism, ".rds")))
  # ## Process genesets from Gene ontology
  # message(format(Sys.time(), " Downloading genesets from Gene Ontology ..."))
  # avail_gaf = c("goa_human.gaf.gz", "mgi.gaf.gz"); names(avail_gaf) = c("hsa", "mmu")
  # base_url = "http://geneontology.org/gene-associations/"
  # locfname = file.path(system.file("extdata", package = "MAGeCKFlute"), "tmp.gaf.gz")
  # download.file(paste0(base_url, avail_gaf[organism]), locfname, quiet = TRUE)
  # go <- read.table(locfname, sep = "\t", comment.char = "!", quote = "", stringsAsFactors = FALSE)
  # file.remove(locfname)
  # colnames(go) <- c("DB","DB_Object_ID","DB_Object_Symbol","Qualifier","GO_ID","DB_Reference(s)",
  #                   "Evidence_Code","With_From","Aspect",
  #                   "DB_Object_Name","DB_Object_Synonym","DB_Object_Type","Taxon",
  #                   "Date","Assigned_By","Annotation_Extension","Gene_Product_Form_ID")
  # go <- go[, c("DB_Object_Symbol", "GO_ID", "Aspect")]
  # idx <- duplicated(paste(go$DB_Object_Symbol, go$GO_ID, go$Aspect, sep = "."))
  # go <- go[!idx, ]
  # go$GO_ID[go$Aspect=="F"] = gsub("GO", "GOMF", go$GO_ID[go$Aspect=="F"])
  # go$GO_ID[go$Aspect=="C"] = gsub("GO", "GOCC", go$GO_ID[go$Aspect=="C"])
  # go$GO_ID[go$Aspect=="P"] = gsub("GO", "GOBP", go$GO_ID[go$Aspect=="P"])
  # tmp = as.data.frame(GO.db::GOTERM)
  # tmp = tmp[!duplicated(tmp$go_id), ]
  # rownames(tmp) = tmp$go_id
  # go$Term = tmp[gsub("GO..", "GO", go$GO_ID), "Term"]
  # go$Entrez = TransGeneID(go$DB_Object_Symbol, "Symbol", "Entrez", organism = organism)
  # go = na.omit(go[, c(5,2,4)])
  # locfname = file.path(system.file("extdata", package = "MAGeCKFlute"),
  #                      paste0("go.all.entrez.", organism, ".rds"))
  # saveRDS(go, locfname)
  # res <- aggregate(go$DB_Object_Symbol, by=list(go$GO_ID), FUN=paste, collapse = "\t")
}
