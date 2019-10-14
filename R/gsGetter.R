#' Extract pathway annotation from GMT file.
#'
#' @docType methods
#' @name gsGetter
#' @rdname gsGetter
#'
#' @param gmtpath The path to customized gmt file.
#' @param type Molecular signatures for testing, available datasets include
#' Pathway (PID, KEGG, REACTOME, BIOCARTA, C2CP), GO (GOBP, GOCC, GOMF),
#' Complex (CORUM, CPX), c1, c2, c3, c4, c6, c7, HALLMARK. It also allows any
#' combination of them (e.g. 'GOBP+GOMF+KEGG+REACTOME').
#' @param limit A two-length vector, specifying the minimal and
#' maximal size of gene sets to load.
#' @param organism 'hsa' or 'mmu'.
#'
#' @return A three-column data frame.
#'
#' @author Wubing Zhang
#'
#' @examples
#' gene2path = gsGetter(type = "REACTOME+CORUM")
#' head(gene2path)
#'
#' @export
#'
gsGetter <- function(gmtpath = NULL, type = "All", limit = c(0, Inf), organism = 'hsa'){
  ## Normalize type
  type = toupper(unlist(strsplit(type, "\\+")))
  if("PATHWAY" %in% type) type = c("PID", "KEGG", "REACTOME", "BIOCARTA", "C2CP", type)
  if("GO" %in% type) type = c("GOBP", "GOCC", "GOMF", type)
  if("COMPLEX" %in% type) type = c("CORUM", "CPX", type)

  ## read GMT files
  if(!is.null(gmtpath)){
    gene2path = ReadGMT(gmtpath, limit = limit)
  }else{
    gene2path = data.frame()
    msigdb = c("hsa_c1.all.v7.0.entrez.gmt.gz", "hsa_c2.cgp.v7.0.entrez.gmt.gz",
               "hsa_c3.all.v7.0.entrez.gmt.gz", "hsa_c4.all.v7.0.entrez.gmt.gz",
               "hsa_c6.all.v7.0.entrez.gmt.gz", "hsa_c7.all.v7.0.entrez.gmt.gz",
               "hsa_go.all.v7.0.entrez.gmt.gz", "hsa_h.all.v7.0.entrez.gmt.gz",
               "hsa_pathway_entrez.gmt.gz", "hsa_complex_entrez.gmt.gz")
    names(msigdb) = toupper(gsub("hsa_|\\..*|_.*", "", msigdb))
    map = c(paste0("C", 1:7), "H", rep("GO",3), rep("PATHWAY", 5), rep("COMPLEX",2))
    names(map) = c(paste0("C", 1:7), "H", "GOBP", "GOCC", "GOMF",
                   "PID", "KEGG", "REACTOME", "BIOCARTA", "C2CP", "CORUM", "CPX")
    if(organism == "hsa"){
      tmp = unique(msigdb[map[type]]); tmp = tmp[!is.na(tmp)]
      msigfiles = file.path(system.file("extdata", package = "MAGeCKFlute"), tmp)
      for(f in msigfiles){
        gene2path = rbind(gene2path, ReadGMT(f, limit = limit))
      }
    }else{
      msigfiles = list.files(system.file("extdata", package = "MAGeCKFlute"),
                             paste0(organism, ".*_entrez.gmt.gz"), full.names = TRUE)
    }
  }
  names(gene2path) = c("Gene","PathwayID", "PathwayName")
  gene2path$PathwayName = toupper(gsub("_", " ", gene2path$PathwayName))
  if(!"All" %in% type){
    idx = toupper(gsub("_.*", "", gene2path$PathwayID)) %in% type
    if(sum(idx)>0) gene2path = gene2path[idx, ]
  }
  gene2path$Gene = as.character(gene2path$Gene)
  gene2path = gene2path[!is.na(gene2path$Gene), ]
  return(gene2path)
}
