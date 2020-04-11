#' Kegg pathway view
#'
#' Plot kegg pathway and color specific genes.
#'
#' @docType methods
#' @name KeggPathwayView
#' @rdname KeggPathwayView
#'
#' @param gene.data Either vector (single sample) or a matrix-like data (multiple sample).
#' Vector should be numeric with gene IDs as names or it may also be character of gene IDs.
#' Character vector is treated as discrete or count data. Matrix-like data structure has
#' genes as rows and samples as columns. Row names should be gene IDs. Here gene ID is a
#' generic concepts, including multiple types of gene, transcript and protein uniquely
#' mappable to KEGG gene IDs. KEGG ortholog IDs are also treated as gene IDs as to handle
#' metagenomic data. Check details for mappable ID types. Default gene.data=NULL.
#' @param cpd.data The same as gene.data, excpet named with IDs mappable to KEGG compound
#' IDs. Over 20 types of IDs included in CHEMBL database can be used here. Check details
#' for mappable ID types. Default cpd.data=NULL. Note that gene.data and cpd.data can't
#' be NULL simultaneously.
#' @param pathway.id Character vector, the KEGG pathway ID(s), usually 5 digit, may also
#' include the 3 letter KEGG species code.
#' @param species Character, either the kegg code, scientific name or the common name of
#' the target species. This applies to both pathway and gene.data or cpd.data. When KEGG
#' ortholog pathway is considered, species="ko". Default species="hsa", it is equivalent
#' to use either "Homo sapiens" (scientific name) or "human" (common name).
#' @param kegg.dir Character, the directory of KEGG pathway data file (.xml) and image file
#'  (.png). Users may supply their own data files in the same format and naming convention
#'   of KEGG's (species code + pathway id, e.g. hsa04110.xml, hsa04110.png etc) in this
#'   directory. Default kegg.dir="." (current working directory).
#' @param cpd.idtype Character, ID type used for the cpd.data. Default cpd.idtype="kegg"
#' (include compound, glycan and drug accessions).
#' @param gene.idtype Character, ID type used for the gene.data, case insensitive. Default
#' gene.idtype="entrez", i.e. Entrez Gene, which are the primary KEGG gene ID for many
#' common model organisms. For other species, gene.idtype should be set to "KEGG" as KEGG
#' use other types of gene IDs. For the common model organisms, you may also specify other
#' types of valid IDs. To check the ID list,
#'  do: data(gene.idtype.list); gene.idtype.list.
#' @param gene.annotpkg Character, the name of the annotation package to use for mapping
#' between other gene ID types including symbols and Entrez gene ID. Default gene.annotpkg=NULL.
#' @param min.nnodes Integer, minimal number of nodes of type "gene","enzyme", "compound"
#' or "ortholog" for a pathway to be considered. Default min.nnodes=3.
#' @param kegg.native Logical, whether to render pathway graph as native KEGG graph (.png)
#'  or using graphviz layout engine (.pdf). Default kegg.native=TRUE.
#' @param map.null Logical, whether to map the NULL gene.data or cpd.data to pathway.
#' When NULL data are mapped, the gene or compound nodes in the pathway will be rendered
#' as actually mapped nodes, except with NA-valued color. When NULL data are not mapped,
#' the nodes are rendered as unmapped nodes. This argument mainly affects native KEGG
#' graph view, i.e. when kegg.native=TRUE. Default map.null=TRUE.
#' @param expand.node Logical, whether the multiple-gene nodes are expanded into
#' single-gene nodes. Each expanded single-gene nodes inherits all edges from the original
#'  multiple-gene node. This option only affects graphviz graph view, i.e. when
#'  kegg.native=FALSE. This option is not effective for most metabolic pathways where it
#'  conflits with converting reactions to edges. Default expand.node=FLASE.
#' @param split.group Logical, whether split node groups are split to individual nodes.
#' Each split member nodes inherits all edges from the node group. This option only affects
#'  graphviz graph view, i.e. when kegg.native=FALSE. This option also effects most
#'  metabolic pathways even without group nodes defined orginally. For these pathways,
#'  genes involved in the same reaction are grouped automatically when converting reactions
#'   to edges unless split.group=TRUE. d split.group=FLASE.
#' @param map.symbol Logical, whether map gene IDs to symbols for gene node labels or use
#' the graphic name from the KGML file. This option is only effective for kegg.native=FALSE
#'  or same.layer=FALSE when kegg.native=TRUE. For same.layer=TRUE when kegg.native=TRUE,
#'  the native KEGG labels will be kept. Default map.symbol=TRUE.
#' @param map.cpdname Logical, whether map compound IDs to formal names for compound node
#' labels or use the graphic name from the KGML file (KEGG compound accessions). This
#' option is only effective for kegg.native=FALSE. When kegg.native=TRUE, the native KEGG
#' labels will be kept. Default map.cpdname=TRUE.
#' @param node.sum Character, the method name to calculate node summary given that multiple
#'  genes or compounds are mapped to it. Poential options include "sum","mean", "median",
#'  "max", "max.abs" and "random". Default node.sum="sum".
#' @param discrete A list of two logical elements with "gene" and "cpd" as the names. This
#'  argument tells whether gene.data or cpd.data should be treated as discrete. Default
#'  dsicrete=list(gene=FALSE, cpd=FALSE), i.e. both data should be treated as continuous.
#' @param limit A list of two numeric elements with "gene" and "cpd" as the names. This
#' argument specifies the limit values for gene.data and cpd.data when converting them to
#' pseudo colors. Each element of the list could be of length 1 or 2. Length 1 suggests
#' discrete data or 1 directional (positive-valued) data, or the absolute limit for 2
#' directional data. Length 2 suggests 2 directional data. Default limit=list(gene=1, cpd=1).
#' @param bins A list of two integer elements with "gene" and "cpd" as the names. This
#' argument specifies the number of levels or bins for gene.data and cpd.data when
#' converting them to pseudo colors. Default limit=list(gene=10, cpd=10).
#' @param both.dirs A list of two logical elements with "gene" and "cpd" as the names.
#' This argument specifies whether gene.data and cpd.data are 1 directional or 2 directional
#' data when converting them to pseudo colors. Default limit=list(gene=TRUE, cpd=TRUE).
#' @param trans.fun A list of two function (not character) elements with "gene" and "cpd"
#' as the names. This argument specifies whether and how gene.data and cpd.data are transformed.
#' Examples are log, abs or users' own functions. Default limit=list(gene=NULL, cpd=NULL).
#' @param low A list of two colors with "gene" and "cpd" as the names.
#' @param mid A list of two colors with "gene" and "cpd" as the names.
#' @param high A list of two colors with "gene" and "cpd" as the names.
#' @param na.col Color used for NA's or missing values in gene.data and cpd.data. d na.col="transparent".
#' @param verbose Boolean
#' @param \dots Extra arguments passed to keggview.native or keggview.graph function.
#'
#' @details The function KeggPathwayView is a revised version of pathview function in pathview package.
#' KeggPathwayView maps and renders user data on relevant pathway graphs. KeggPathwayView
#' is a stand alone program for pathway based data integration and visualization. It also
#' seamlessly integrates with pathway and functional analysis tools for large-scale and
#' fully automated analysis. KeggPathwayView provides strong support for data Integration.
#' It works with: 1) essentially all types of biological data mappable to pathways, 2) over
#' 10 types of gene or protein IDs, and 20 types of compound or metabolite IDs, 3) pathways
#' for over 2000 species as well as KEGG orthology, 4) varoius data attributes and formats,
#' i.e. continuous/discrete data, matrices/vectors, single/multiple samples etc.
#' To see mappable external gene/protein IDs do: data(gene.idtype.list), to see mappable
#' external compound related IDs do: data(rn.list); names(rn.list). KeggPathwayView
#' generates both native KEGG view and Graphviz views for pathways. Currently only KEGG
#' pathways are implemented. Hopefully, pathways from Reactome, NCI and other databases
#' will be supported in the future.
#'
#' The argument \code{low}, \code{mid}, and \code{high} specifies the color spectra to code gene.data and cpd.data. When data are
#'  1 directional (TRUE value in both.dirs), only mid and high are used to specify the
#'  color spectra. Default spectra (low-mid-high) "green"-"gray"-"red" and "blue"-"gray"-"yellow"
#'  are used for gene.data and cpd.data respectively. The values for 'low, mid, high' can
#'  be given as color names ('red'), plot color index (2=red), and HTML-style RGB, ("\#FF0000"=red).

#'
#' @return The result returned by KeggPathwayView function is a named list corresponding to
#' the input pathway ids. Each element (for each pathway itself is a named list, with 2 elements
#' ("plot.data.gene", "plot.data.cpd"). Both elements are data.frame or NULL depends on the
#' corresponding input data gene.data and cpd.data. These data.frames record the plot data for
#' mapped gene or compound nodes: rows are mapped genes/compounds, columns are:
#' \item{kegg.names }{standard KEGG IDs/Names for mapped nodes. It's Entrez Gene ID or KEGG
#' Compound Accessions.}
#' \item{labels }{Node labels to be used when needed.}
#' \item{all.mapped }{All molecule (gene or compound) IDs mapped to this node.}
#' \item{type }{node type, currently 4 types are supported: "gene","enzyme", "compound" and
#' "ortholog".}
#' \item{x }{x coordinate in the original KEGG pathway graph.}
#' \item{y }{y coordinate in the original KEGG pathway graph.}
#' \item{width }{node width in the original KEGG pathway graph.}
#' \item{height }{node height in the original KEGG pathway graph.}
#' \item{other columns}{columns of the mapped gene/compound data and corresponding
#' pseudo-color codes for individual samples}
#'
#' @author Wubing Zhang
#'
#'
#' @examples
#' #load data
#' data(mle.gene_summary)
#' dd = ReadBeta(mle.gene_summary)
#' gene.data = dd$plx
#' names(gene.data) = rownames(dd)
#'
#' pv.out <- KeggPathwayView(gene.data, pathway.id = "04110",
#'   species = "hsa", out.suffix = "gse16873", kegg.native = TRUE)
#'
#' @importFrom KEGGREST keggConv
#' @import pathview
#' @export
KeggPathwayView=function (gene.data = NULL, cpd.data = NULL, pathway.id,
                          species = "hsa", kegg.dir = ".",
                          cpd.idtype = "kegg",gene.idtype ="ENTREZ",
                          gene.annotpkg = NULL, min.nnodes = 3,
                          kegg.native = TRUE, map.null = TRUE,
                          expand.node = FALSE, split.group = FALSE,
                          map.symbol = TRUE, map.cpdname = TRUE,
                          node.sum = "sum", discrete =
                            list(gene = FALSE, cpd = FALSE),
                          limit=list(gene=1, cpd=1),
                          bins = list(gene = 10, cpd = 10),
                          both.dirs = list(gene =TRUE,cpd =TRUE),
                          trans.fun = list(gene = NULL, cpd = NULL),
                          low = list(gene = "deepskyblue1", cpd = "blue"),
                          mid = list(gene = "gray", cpd = "gray"),
                          high = list(gene = "red", cpd ="yellow"),
                          na.col = "transparent", verbose = TRUE, ...)
{
  requireNamespace("KEGGREST")
  dtypes = !is.null(gene.data) + (!is.null(cpd.data))
  cond0 = dtypes == 1 & is.numeric(limit) & length(limit) > 1
  if (cond0) {
    if (limit[1] != limit[2] & is.null(names(limit)))
      limit = list(gene = limit[1:2], cpd = limit[1:2])
  }
  if (is.null(trans.fun)) trans.fun = list(gene = NULL, cpd = NULL)
  arg.len2 = c("discrete", "limit", "bins", "both.dirs", "trans.fun",
               "low", "mid", "high")
  for (arg in arg.len2) {
    obj1 = eval(as.name(arg))
    if (length(obj1) == 1)
      obj1 = rep(obj1, 2)
    if (length(obj1) > 2)
      obj1 = obj1[1:2]
    obj1 = as.list(obj1)
    ns = names(obj1)
    if (length(ns) == 0 | !all(c("gene", "cpd") %in% ns))
      names(obj1) = c("gene", "cpd")
    assign(arg, obj1)
  }
  if (is.character(gene.data)) {
    gd.names = gene.data
    gene.data = rep(1, length(gene.data))
    names(gene.data) = gd.names
    both.dirs$gene = FALSE
    ng = length(gene.data)
    nsamp.g = 1
  }else if (!is.null(gene.data)) {
    if (length(dim(gene.data)) == 2) {
      gd.names = rownames(gene.data)
      ng = nrow(gene.data)
      nsamp.g = 2
    }else if (is.numeric(gene.data) & is.null(dim(gene.data))) {
      gd.names = names(gene.data)
      ng = length(gene.data)
      nsamp.g = 1
    }else stop("wrong gene.data format!")
  }else if (is.null(cpd.data)) {
    stop("gene.data and cpd.data are both NULL!")
  }
  gene.idtype = toupper(gene.idtype)
  bods <- data.frame(package = paste0("org.", c("Hs", "Mm", "Rn", "Bt", "Cf", "Pt", "Ss"), ".eg.db"),
                     species = c("Human", "Mouse", "Rat", "Bovine", "Canine", "Chimp", "Pig"),
                     "kegg code" = c("hsa", "mmu", "rno", "bta", "cfa", "ptr", "ssc"),
                     check.names = FALSE, stringsAsFactors = FALSE)
  # if (species != "ko") {
  #   species.data = kegg.species.code(species, na.rm =TRUE,
  #                                    code.only = FALSE)
  # }else {
  #   species.data = c(kegg.code = "ko", entrez.gnodes = "0",
  #                    kegg.geneid = "K01488", ncbi.geneid = "")
  #   gene.idtype = "KEGG"
  #   msg.fmt = "Only KEGG ortholog gene ID is supported,
  #     make sure it looks like \"%s\"!"
  #   msg = sprintf(msg.fmt, species.data["kegg.geneid"])
  #   if(verbose) message(Sys.time(), " # Note: ", msg)
  # }
  # if (length(dim(species.data)) == 2) {
  #   if(verbose) message(Sys.time(), " # Note: More than two valide species!")
  #   species.data = species.data[1, ]
  # }
  # species = species.data["kegg.code"]
  # entrez.gnodes = species.data["entrez.gnodes"] == 1
  entrez.gnodes = TRUE
  # if (is.na(species.data["ncbi.geneid"])) {
  #   if (!is.na(species.data["kegg.geneid"])) {
  #     msg.fmt = "Only native KEGG gene ID is supported for this species,
  #     \nmake sure it looks like \"%s\"!"
  #     msg = sprintf(msg.fmt, species.data["kegg.geneid"])
  #     if(verbose) message(Sys.time(), " #  Note: ", msg)
  #   }else {
  #     stop("This species is not annotated in KEGG!")
  #   }
  # }
  if (is.null(gene.annotpkg))
    gene.annotpkg = bods[match(species, bods[, 3]), 1]
  if (length(grep("ENTREZ|KEGG", gene.idtype)) < 1 & !is.null(gene.data)) {
    if (is.na(gene.annotpkg))
      stop("No proper gene annotation package available!")
    if (!gene.idtype %in% c("ENTREZ|SYMBOL|KEGG"))
      stop("Wrong input gene ID type!")
    gene.idmap = TransGeneID(gd.names, fromType = gene.idtype,
                             toType = "ENTREZ", organism = species)
    gene.idmap = cbind(names(gene.idmap), gene.idmap)
    gene.data = mol.sum(gene.data, gene.idmap)
    gene.idtype = "ENTREZ"
  }
  if (gene.idtype == "ENTREZ" & !entrez.gnodes & !is.null(gene.data)) {
    if(verbose) message(Sys.time(), " #  Info: Getting gene ID data from KEGG...")
    gene.idmap = keggConv("ncbi-geneid", species)
    if(verbose) message(Sys.time(), " #  Info: Done with data retrieval!")
    kegg.ids = gsub(paste(species, ":", sep = ""), "", names(gene.idmap))
    ncbi.ids = gsub("ncbi-geneid:", "", gene.idmap)
    gene.idmap = cbind(ncbi.ids, kegg.ids)
    gene.data = mol.sum(gene.data, gene.idmap)
    gene.idtype = "KEGG"
  }
  if (is.character(cpd.data)) {
    cpdd.names = cpd.data
    cpd.data = rep(1, length(cpd.data))
    names(cpd.data) = cpdd.names
    both.dirs$cpd = FALSE
    ncpd = length(cpd.data)
  }else if (!is.null(cpd.data)) {
    if (length(dim(cpd.data)) == 2) {
      cpdd.names = rownames(cpd.data)
      ncpd = nrow(cpd.data)
    }else if (is.numeric(cpd.data) & is.null(dim(cpd.data))) {
      cpdd.names = names(cpd.data)
      ncpd = length(cpd.data)
    }else stop("wrong cpd.data format!")
  }
  # if (length(grep("kegg", cpd.idtype)) < 1 & !is.null(cpd.data)) {
  #   rn.list = NULL
  #   data(rn.list, package = "pathview")
  #   cpd.types = c(names(rn.list), "name")
  #   cpd.types = tolower(cpd.types)
  #   cpd.types = cpd.types[-grep("kegg", cpd.types)]
  #   if (!tolower(cpd.idtype) %in% cpd.types)
  #     stop("Wrong input cpd ID type!")
  #   cpd.idmap = cpd2kegg(cpdd.names, in.type = cpd.idtype)
  #   cpd.data = mol.sum(cpd.data, cpd.idmap)
  # }
  warn.fmt = "    Parsing %s file failed, please check the file!"
  if (length(grep(species, pathway.id)) > 0) {
    pathway.name = pathway.id
    pathway.id = gsub(species, "", pathway.id)
  }else pathway.name = paste(species, pathway.id, sep = "")
  kfiles = list.files(path = kegg.dir, pattern = "[.]xml|[.]png")
  npath = length(pathway.id)
  out.list = list()
  tfiles.xml = paste(pathway.name, "xml", sep = ".")
  tfiles.png = paste(pathway.name, "png", sep = ".")
  if (kegg.native){
    ttype = c("xml", "png")
  }else ttype = "xml"
  xml.file <- paste(kegg.dir, "/", tfiles.xml, sep = "")
  for (i in 1:npath) {
    if (kegg.native){
      tfiles = c(tfiles.xml[i], tfiles.png[i])
    }else tfiles = tfiles.xml[i]
    if (!all(tfiles %in% kfiles)) {
      dstatus = download.kegg(pathway.id = pathway.id[i],
                              species = species, kegg.dir = kegg.dir,
                              file.type = ttype)
      if (dstatus == "failed") {
        warn.fmt = "    Failed to download KEGG xml/png files, %s skipped!"
        warn.msg = sprintf(warn.fmt, pathway.name[i])
        warning(warn.msg)
        return(invisible(0))
      }
    }
    if (kegg.native) {
      node.data = try(node.info(xml.file[i]), silent =TRUE)
      if (is(node.data, "try-error")) {
        warn.msg = sprintf(warn.fmt, xml.file[i])
        warning(warn.msg)
        return(invisible(0))
      }
      node.type = c("gene", "enzyme", "compound", "ortholog")
      sel.idx = node.data$type %in% node.type
      nna.idx = !is.na(node.data$x + node.data$y + node.data$width +
                         node.data$height)
      sel.idx = sel.idx & nna.idx
      if (sum(sel.idx) < min.nnodes) {
        warn.fmt = "    Number of mappable nodes is below %d, %s skipped!"
        warn.msg = sprintf(warn.fmt, min.nnodes, pathway.name[i])
        warning(warn.msg)
        return(invisible(0))
      }
      node.data = lapply(node.data, "[", sel.idx)
    }else {
      gR1 = try(parseKGML2Graph2(xml.file[i], genes =FALSE,
                                 expand = expand.node, split.group =
                                   split.group),
                silent =TRUE)
      node.data = try(node.info(gR1), silent =TRUE)
      if (is(node.data, "try-error")) {
        warn.msg = sprintf(warn.fmt, xml.file[i])
        warning(warn.msg)
        return(invisible(0))
      }
    }
    if (species == "ko"){
      gene.node.type = "ortholog"
    }else gene.node.type = "gene"
    head(gene.data)
    if ((!is.null(gene.data) | map.null) & sum(node.data$type ==
                                               gene.node.type) > 1) {
      plot.data.gene = node.map(gene.data, node.data,
                                node.types = gene.node.type,
                                node.sum = node.sum, entrez.gnodes = entrez.gnodes)
      plot.data.gene<-plot.data.gene[
        rowSums(plot.data.gene[,c("x","y","width","height")])!=4,]
      kng = plot.data.gene$kegg.names
      kng.char = gsub("[0-9]", "", unlist(kng))
      if (any(kng.char > ""))
        entrez.gnodes = FALSE
      if (map.symbol & species != "ko" & entrez.gnodes) {
        if (is.na(gene.annotpkg)) {
          warn.fmt = "    No annotation package for the species %s,
          gene symbols not mapped!"
          warn.msg = sprintf(warn.fmt, species)
          warning(warn.msg)
        }else {
          #=====My revised===========
          plot.data.gene$labels = TransGeneID(as.character(
            plot.data.gene$kegg.names),"Entrez", "Symbol", organism = species)[as.character(plot.data.gene$kegg.names)]
          #==========================
          mapped.gnodes = rownames(plot.data.gene)
          node.data$labels[mapped.gnodes] = plot.data.gene$labels
        }
      }
      cols.ts.gene = node.color(plot.data.gene, limit$gene,
                                bins$gene, both.dirs = both.dirs$gene,
                                trans.fun = trans.fun$gene,
                                discrete = discrete$gene, low = low$gene,
                                mid = mid$gene,
                                high = high$gene, na.col = na.col)
    }else plot.data.gene = cols.ts.gene = NULL
    if ((!is.null(cpd.data) | map.null) & sum(node.data$type ==
                                              "compound") > 1) {
      plot.data.cpd = node.map(cpd.data, node.data,
                               node.types = "compound",
                               node.sum = node.sum)
      if (map.cpdname & !kegg.native) {
        plot.data.cpd$labels = cpdkegg2name(plot.data.cpd$labels)[,2]
        mapped.cnodes = rownames(plot.data.cpd)
        node.data$labels[mapped.cnodes] = plot.data.cpd$labels
      }
      cols.ts.cpd = node.color(plot.data.cpd, limit$cpd,
                               bins$cpd, both.dirs = both.dirs$cpd,
                               trans.fun = trans.fun$cpd,
                               discrete = discrete$cpd, low = low$cpd,
                               mid = mid$cpd,
                               high = high$cpd, na.col = na.col)
    }else plot.data.cpd = cols.ts.cpd = NULL
    if (kegg.native) {
      pv.pars = keggview.native(plot.data.gene = plot.data.gene,
                                cols.ts.gene = cols.ts.gene,
                                plot.data.cpd = plot.data.cpd,
                                cols.ts.cpd = cols.ts.cpd,
                                node.data = node.data,
                                pathway.name = pathway.name[i],
                                kegg.dir = kegg.dir,
                                limit = limit, bins = bins,
                                both.dirs = both.dirs, discrete = discrete,
                                low = low, mid = mid, high = high,
                                na.col = na.col, ...)
    }else {
      pv.pars = keggview.graph(plot.data.gene = plot.data.gene,
                               cols.ts.gene = cols.ts.gene,
                               plot.data.cpd = plot.data.cpd,
                               cols.ts.cpd = cols.ts.cpd,
                               node.data = node.data,
                               path.graph = gR1,
                               pathway.name = pathway.name[i],
                               map.cpdname = map.cpdname,
                               split.group = split.group,
                               limit = limit, bins = bins,
                               both.dirs = both.dirs, discrete = discrete,
                               low = low, mid = mid, high = high,
                               na.col = na.col)
    }
    plot.data.gene = cbind(plot.data.gene, cols.ts.gene)
    if (!is.null(plot.data.gene)) {
      cnames = colnames(plot.data.gene)[-(1:8)]
      nsamp = length(cnames)/2
      if (nsamp > 1) {
        cnames[(nsamp + 1):(2 * nsamp)] =
          paste(cnames[(nsamp + 1):(2 * nsamp)], "col", sep = ".")
      }
      else cnames[2] = "mol.col"
      colnames(plot.data.gene)[-(1:8)] = cnames
    }
    plot.data.cpd = cbind(plot.data.cpd, cols.ts.cpd)
    if (!is.null(plot.data.cpd)) {
      cnames = colnames(plot.data.cpd)[-(1:8)]
      nsamp = length(cnames)/2
      if (nsamp > 1) {
        cnames[(nsamp + 1):(2 * nsamp)] =
          paste(cnames[(nsamp + 1):(2 * nsamp)], "col", sep = ".")
      }
      else cnames[2] = "mol.col"
      colnames(plot.data.cpd)[-(1:8)] = cnames
    }
    out.list[[i]] = list(plot.data.gene = plot.data.gene,
                         plot.data.cpd = plot.data.cpd)
  }
  if (npath == 1)
    out.list = out.list[[1]]
  else names(out.list) = pathway.name
  return(invisible(out.list))
}

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
#' @param view_allpath boolean, specifying whether view all pathways. Default view_allpath='FALSE',
#'  and only plot top enriched pathways.
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
#' @seealso \code{\link{KeggPathwayView}}
#'
#' @examples
#' data(mle.gene_summary)
#' # Read beta score from gene summary table in MAGeCK MLE results
#' dd = ReadBeta(mle.gene_summary)
#' colnames(dd)[2:3] = c("Control", "Treatment")
#' arrangePathview(dd, "hsa00534", title=NULL, sub=NULL, organism="hsa")
#'
#' @importFrom png readPNG
#' @importFrom grid rasterGrob
#' @importFrom gridExtra grid.arrange
#' @export

arrangePathview <- function(genelist, pathways=c(), top = 4, ncol = 2,
                            title=NULL, sub=NULL,
                            organism='hsa', view_allpath= FALSE,
                            output=".", path.archive = ".",
                            kegg.native = TRUE, verbose = TRUE){
  #====No pathways supplied======================
  if(length(pathways)<1){
    p=ggplot()
    p=p+geom_text(aes(x=0,y=0,label="No enriched terms"),size=6)
    p=p+theme_void()
    print(p)
    return(0)
  }

  if(length(pathways)<top || view_allpath){
    keggID=pathways
  }else{
    keggID=pathways[1:top]
  }

  if(verbose) message(Sys.time(), " # Starting plot kegg pathways for ", sub, " ", title)

  p1 <- suppressWarnings(KeggPathwayView(gene.data  = genelist[,c("Control","Treatment")],
                        pathway.id = keggID, species=organism, kegg.dir = path.archive,
                        kegg.native = kegg.native))

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



