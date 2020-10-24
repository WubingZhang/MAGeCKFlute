#' Kegg pathway view and arrange grobs on page
#'
#' Kegg pathway view and arrange grobs on page.
#'
#' @docType methods
#' @name arrangePathview
#' @rdname arrangePathview
#'
#' @param genelist a data frame with columns of ENTREZID, Control and Treatment. The columns
#'  of Control and Treatment represent gene score in Control and Treatment sample.
#' @param pathways character vector, the KEGG pathway ID(s), usually 5 digit, may also
#' include the 3 letter KEGG species code.
#' @param top integer, specifying how many top enriched pathways to be visualized.
#' @param ncol integer, specifying how many column of figures to be arranged in each page.
#' @param title  optional string, or grob.
#' @param sub optional string, or grob.
#' @param organism character, either the kegg code, scientific name or the common name of
#' the target species. This applies to both pathway and gene.data or cpd.data. When KEGG
#' ortholog pathway is considered, species="ko". Default species="hsa", it is equivalent
#' to use either "Homo sapiens" (scientific name) or "human" (common name).
#' @param output Path to save plot to.
#' @param path.archive character, the directory of KEGG pathway data file (.xml) and image file
#'  (.png). Users may supply their own data files in the same format and naming convention
#'   of KEGG's (species code + pathway id, e.g. hsa04110.xml, hsa04110.png etc) in this
#'   directory. Default kegg.dir="." (current working directory).
#' @param kegg.native logical, whether to render pathway graph as native KEGG graph (.png)
#'  or using graphviz layout engine (.pdf). Default kegg.native=TRUE.
#' @param verbose Boolean
#'
#' @return plot on the current device
#'
#' @author Wubing Zhang
#'
#' @examples
#' file3 = file.path(system.file("extdata", package = "MAGeCKFlute"),
#' "testdata/mle.gene_summary.txt")
#' dd = ReadBeta(file3)
#' colnames(dd)[2:3] = c("Control", "Treatment")
#' # arrangePathview(dd, c("hsa00534"), title=NULL, sub=NULL, organism="hsa")
#'
#' @importFrom grid rasterGrob
#' @importFrom gridExtra grid.arrange
#' @export

arrangePathview <- function(genelist, pathways=c(), top = 4, ncol = 2,
                            title=NULL, sub=NULL,
                            organism='hsa',
                            output=".", path.archive = ".",
                            kegg.native = TRUE, verbose = TRUE){
  if (!requireNamespace("png", quietly = TRUE)) {
    stop("Package \"png\" is required. Please install it.", call. = FALSE)
  }
  if (!require("pathview", quietly = TRUE)) {
    stop("Package \"pathview\" is required. Please install it.", call. = FALSE)
  }
  #====No pathways supplied======================
  if(length(pathways)<1){
    p=ggplot()
    p=p+geom_text(aes(x=0,y=0,label="No enriched terms"),size=6)
    p=p+theme_void()
    print(p)
    return(0)
  }

  if(length(pathways)<=top){
    keggID=pathways
  }else{
    keggID=pathways[1:top]
  }
  if(length(keggID)<1) return()

  if(verbose) message(Sys.time(), " # Starting plot kegg pathways for ", sub, " ", title)
  for(id in keggID){
    tryCatch(pathview(gene.data  = genelist[,c("Control","Treatment")],
                 pathway.id = id, species=organism, kegg.dir = path.archive,
                 kegg.native = kegg.native), error = function(e)
                   message("pathview failed for ", id))
  }
  #Maybe there are not multi file, but only keggID.pathview.png
  allpngnames=paste0(keggID, ".pathview.multi.png")
  idx = file.exists(allpngnames)
  allpngnames=allpngnames[idx]

  if(length(allpngnames)>0){
    toFile=file.path(output, paste0(title,"_",sub,"_",allpngnames))
    boo=file.rename(from=allpngnames,to=toFile)
  }else{boo=FALSE}
  originPng=paste0(keggID, ".png")
  originXML=paste0(keggID, ".xml")
  failMulti=paste0(keggID, ".pathview.png")
  suppressWarnings(file.remove(originPng))
  suppressWarnings(file.remove(originXML))
  suppressWarnings(file.remove(failMulti))

  if(all(boo)){
    pngnames = file.path(output, paste0(title, "_", sub, "_", allpngnames))
    idx=file.exists(pngnames)
    pngnames = pngnames[idx]
  }else{pngnames=c()}

  if(length(pngnames)>0){
    thePlots <- lapply (pngnames, function(figure) {
      grid::rasterGrob(png::readPNG(figure, native = FALSE),interpolate = FALSE)})
    for(i in 1:ceiling(length(thePlots)/ncol)){
      if(ncol*i <= length(thePlots))
        do.call(gridExtra::grid.arrange, c(thePlots[(ncol*(i-1)+1):(ncol*i)], ncol = ncol, top=title, bottom=sub))
      else
        do.call(gridExtra::grid.arrange, c(thePlots[(ncol*(i-1)+1):length(thePlots)], ncol = ncol, top=title, bottom=sub))
    }
  }
}



