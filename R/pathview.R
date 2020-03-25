parseKGML2Graph2 <-function (file, ...)
{
  pathway <- parseKGML2(file)
  gR <- KEGGpathway2Graph2(pathway, ...)
  return(gR)
}

#' @importFrom XML xmlTreeParse xmlRoot xmlChildren xmlName
#' @importFrom KEGGgraph parsePathwayInfo parseEntry
parseKGML2<-function (file)
{
  doc <- XML::xmlTreeParse(file, getDTD = FALSE)
  r <- XML::xmlRoot(doc)
  childnames <- sapply(XML::xmlChildren(r), XML::xmlName)
  isEntry <- childnames == "entry"
  isRelation <- childnames == "relation"
  isReaction <- childnames == "reaction"
  kegg.pathwayinfo <- KEGGgraph::parsePathwayInfo(r)
  kegg.nodes <- sapply(r[isEntry], KEGGgraph::parseEntry)
  kegg.edges <- sapply(r[isRelation], KEGGgraph::parseRelation)
  kegg.reactions <- sapply(r[isReaction], parseReaction2)
  names(kegg.nodes) <- sapply(kegg.nodes, KEGGgraph::getEntryID)
  pathway <- new("KEGGPathway", pathwayInfo = kegg.pathwayinfo,
                 nodes = kegg.nodes, edges = kegg.edges, reactions = kegg.reactions)
  return(pathway)
}
#' @importFrom KEGGgraph splitKEGGgroup expandKEGGPathway nodes edges getEntryID getEntryID subGraphByNodeType
#' @importFrom graph nodeDataDefaults edgeDataDefaults
KEGGpathway2Graph2 <-
  function (pathway, genesOnly = TRUE, expandGenes = TRUE, split.group=FALSE, check.reaction=TRUE)
  {
    requireNamespace("pathview")
    requireNamespace("KEGGgraph")
    stopifnot(is(pathway, "KEGGPathway"))
    if(split.group) pathway <- KEGGgraph::splitKEGGgroup(pathway)
    rdata=(pathway@reactions)
    if (expandGenes){
      if(check.reaction & length(rdata)>0) message("Note: ", "Gene nodes not expanded when reactions are converted to edges!")
      else pathway <- KEGGgraph::expandKEGGPathway(pathway)
    }
    knodes <- KEGGgraph::nodes(pathway)
    kedges <- KEGGgraph::edges(pathway)
    node.entryIDs <- KEGGgraph::getEntryID(knodes)
    edge.entryIDs <- KEGGgraph::getEntryID(kedges)
    V <- node.entryIDs
    edL <- vector("list", length = length(V))
    names(edL) <- V
    if (is.null(nrow(edge.entryIDs))) {
      for (i in seq(along = edL)) {
        edL[[i]] <- list()
      }
    }
    else {
      for (i in 1:length(V)) {
        id <- node.entryIDs[i]
        hasRelation <- id == edge.entryIDs[, "Entry1ID"]
        if (!any(hasRelation)) {
          edL[[i]] <- list(edges = NULL)
        }
        else {
          entry2 <- unname(unique(edge.entryIDs[hasRelation,
                                                "Entry2ID"]))
          edL[[i]] <- list(edges = entry2)
        }
      }
    }
    gR <- new("graphNEL", nodes = V, edgeL = edL, edgemode = "directed")

    if(check.reaction & length(rdata)>0){
      r2e.res=reaction2edge(pathway, gR)
      gR=r2e.res[[1]]
      kedges=r2e.res[[2]]
      knodes=r2e.res[[3]]
    }

    names(kedges) <- sapply(kedges, function(x) paste(getEntryID(x),
                                                      collapse = "~"))
    env.node <- new.env()
    env.edge <- new.env()
    assign("nodes", knodes, envir = env.node)
    assign("edges", kedges, envir = env.edge)
    graph::nodeDataDefaults(gR, "KEGGNode") <- env.node
    graph::edgeDataDefaults(gR, "KEGGEdge") <- env.edge
    if (genesOnly) {
      gR <- KEGGgraph::subGraphByNodeType(gR, "gene")
    }
    return(gR)
  }

#' @importFrom XML xmlAttrs xmlChildren
#'
parseReaction2 <-
  function (reaction)
  {
    attrs <- XML::xmlAttrs(reaction)
    name <- attrs[["name"]]
    type <- attrs[["type"]]
    children <- XML::xmlChildren(reaction)
    childrenNames <- names(children)
    substrateIndices <- grep("^substrate$", childrenNames)
    productIndices <- grep("^product$", childrenNames)
    substrateName <- substrateAltName <- vector("character",
                                                length(substrateIndices))
    productName <- productAltName <- vector("character", length(productIndices))
    for (i in seq(along = substrateIndices)) {
      ind <- substrateIndices[i]
      substrate <- children[[ind]]
      substrateName[i] <- XML::xmlAttrs(substrate)[["id"]]
      substrateChildren <- XML::xmlChildren(substrate)
      if (length(substrateChildren) > 0) {
        substrateAlt <- substrateChildren$alt
        substrateAltName[i] <- XML::xmlAttrs(substrateAlt)[["name"]]
      }
      else {
        substrateAlt <- as.character(NA)
        substrateAltName[i] <- as.character(NA)
      }
    }
    for (i in seq(along = productIndices)) {
      ind <- productIndices[i]
      product <- children[[ind]]
      productName[i] <- XML::xmlAttrs(product)[["id"]]
      productChildren <- XML::xmlChildren(product)
      if (length(productChildren) > 0) {
        productAlt <- productChildren$alt
        productAltName[i] <- XML::xmlAttrs(productAlt)[["name"]]
      }
      else {
        productAlt <- as.character(NA)
        productAltName[i] <- as.character(NA)
      }
    }
    new("KEGGReaction", name = name, type = type, substrateName = substrateName,
        substrateAltName = substrateAltName, productName = productName,
        productAltName = productAltName)
  }
