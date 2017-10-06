\name{KeggPathwayView}
\alias{KeggPathwayView}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot kegg pathway
}
\description{
kegg pathway view
}
\usage{
KeggPathwayView(gene.data = NULL, cpd.data = NULL, pathway.id,
                species = "hsa", kegg.dir = ".", cpd.idtype = "kegg",
                gene.idtype = "ENTREZ", gene.annotpkg = NULL,
                min.nnodes = 3, kegg.native = TRUE, map.null = TRUE,
                expand.node = FALSE, split.group = FALSE,
                map.symbol = TRUE, map.cpdname = TRUE, node.sum = "sum",
                discrete = list(gene = FALSE, cpd = FALSE),
                limit = list(gene = 1, cpd = 1),
                bins = list(gene = 10, cpd = 10),
                both.dirs = list(gene = T, cpd = T),
                trans.fun = list(gene = NULL, cpd = NULL),
                low = list(gene = "deepskyblue1", cpd = "blue"),
                mid = list(gene = "gray", cpd = "gray"),
                high = list(gene = "red", cpd = "yellow"),
                na.col = "transparent", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene.data}{
either vector (single sample) or a matrix-like data (multiple sample). Vector should be numeric with gene IDs as names or it may also be character of gene IDs. Character vector is treated as discrete or count data. Matrix-like data structure has genes as rows and samples as columns. Row names should be gene IDs. Here gene ID is a generic concepts, including multiple types of gene, transcript and protein uniquely mappable to KEGG gene IDs. KEGG ortholog IDs are also treated as gene IDs as to handle metagenomic data. Check details for mappable ID types. Default gene.data=NULL.

}
  \item{cpd.data}{
the same as gene.data, excpet named with IDs mappable to KEGG compound IDs. Over 20 types of IDs included in CHEMBL database can be used here. Check details for mappable ID types. Default cpd.data=NULL. Note that gene.data and cpd.data can't be NULL simultaneously.
}
  \item{pathway.id}{
character vector, the KEGG pathway ID(s), usually 5 digit, may also include the 3 letter KEGG species code.
}
  \item{species}{
character, either the kegg code, scientific name or the common name of the target species. This applies to both pathway and gene.data or cpd.data. When KEGG ortholog pathway is considered, species="ko". Default species="hsa", it is equivalent to use either "Homo sapiens" (scientific name) or "human" (common name).
}
  \item{kegg.dir}{
character, the directory of KEGG pathway data file (.xml) and image file (.png). Users may supply their own data files in the same format and naming convention of KEGG's (species code + pathway id, e.g. hsa04110.xml, hsa04110.png etc) in this directory. Default kegg.dir="." (current working directory).
}
  \item{cpd.idtype}{
character, ID type used for the cpd.data. Default cpd.idtype="kegg" (include compound, glycan and drug accessions).
}
  \item{gene.idtype}{
character, ID type used for the gene.data, case insensitive. Default gene.idtype="entrez", i.e. Entrez Gene, which are the primary KEGG gene ID for many common model organisms. For other species, gene.idtype should be set to "KEGG" as KEGG use other types of gene IDs. For the common model organisms (to check the list, do: data(bods); bods), you may also specify other types of valid IDs. To check the ID list, do: data(gene.idtype.list); gene.idtype.list.
}
  \item{gene.annotpkg}{
character, the name of the annotation package to use for mapping between other gene ID types including symbols and Entrez gene ID. Default gene.annotpkg=NULL.
}
  \item{min.nnodes}{
integer, minimal number of nodes of type "gene","enzyme", "compound" or "ortholog" for a pathway to be considered. Default min.nnodes=3.
}
  \item{kegg.native}{
logical, whether to render pathway graph as native KEGG graph (.png) or using graphviz layout engine (.pdf). Default kegg.native=TRUE.
}
  \item{map.null}{
logical, whether to map the NULL gene.data or cpd.data to pathway. When NULL data are mapped, the gene or compound nodes in the pathway will be rendered as actually mapped nodes, except with NA-valued color. When NULL data are not mapped, the nodes are rendered as unmapped nodes. This argument mainly affects native KEGG graph view, i.e. when kegg.native=TRUE. Default map.null=TRUE.
}
  \item{expand.node}{
logical, whether the multiple-gene nodes are expanded into single-gene nodes. Each expanded single-gene nodes inherits all edges from the original multiple-gene node. This option only affects graphviz graph view, i.e. when kegg.native=FALSE. This option is not effective for most metabolic pathways where it conflits with converting reactions to edges. Default expand.node=FLASE.
}
  \item{split.group}{
logical, whether split node groups are split to individual nodes. Each split member nodes inherits all edges from the node group. This option only affects graphviz graph view, i.e. when kegg.native=FALSE. This option also effects most metabolic pathways even without group nodes defined orginally. For these pathways, genes involved in the same reaction are grouped automatically when converting reactions to edges unless split.group=TRUE. d split.group=FLASE.
}
  \item{map.symbol}{
logical, whether map gene IDs to symbols for gene node labels or use the graphic name from the KGML file. This option is only effective for kegg.native=FALSE or same.layer=FALSE when kegg.native=TRUE. For same.layer=TRUE when kegg.native=TRUE, the native KEGG labels will be kept. Default map.symbol=TRUE.
}
  \item{map.cpdname}{
logical, whether map compound IDs to formal names for compound node labels or use the graphic name from the KGML file (KEGG compound accessions). This option is only effective for kegg.native=FALSE. When kegg.native=TRUE, the native KEGG labels will be kept. Default map.cpdname=TRUE.
}
  \item{node.sum}{
character, the method name to calculate node summary given that multiple genes or compounds are mapped to it. Poential options include "sum","mean", "median", "max", "max.abs" and "random". Default node.sum="sum".
}
  \item{discrete}{
a list of two logical elements with "gene" and "cpd" as the names. This argument tells whether gene.data or cpd.data should be treated as discrete. Default dsicrete=list(gene=FALSE, cpd=FALSE), i.e. both data should be treated as continuous.
}
  \item{limit}{
a list of two numeric elements with "gene" and "cpd" as the names. This argument specifies the limit values for gene.data and cpd.data when converting them to pseudo colors. Each element of the list could be of length 1 or 2. Length 1 suggests discrete data or 1 directional (positive-valued) data, or the absolute limit for 2 directional data. Length 2 suggests 2 directional data. Default limit=list(gene=1, cpd=1).
}
  \item{bins}{
a list of two integer elements with "gene" and "cpd" as the names. This argument specifies the number of levels or bins for gene.data and cpd.data when converting them to pseudo colors. Default limit=list(gene=10, cpd=10).
}
  \item{both.dirs}{
a list of two logical elements with "gene" and "cpd" as the names. This argument specifies whether gene.data and cpd.data are 1 directional or 2 directional data when converting them to pseudo colors. Default limit=list(gene=TRUE, cpd=TRUE).
}
  \item{trans.fun}{
a list of two function (not character) elements with "gene" and "cpd" as the names. This argument specifies whether and how gene.data and cpd.data are transformed. Examples are log, abs or users' own functions. Default limit=list(gene=NULL, cpd=NULL).
}
  \item{low, mid, high}{
each is a list of two colors with "gene" and "cpd" as the names. This argument specifies the color spectra to code gene.data and cpd.data. When data are 1 directional (TRUE value in both.dirs), only mid and high are used to specify the color spectra. Default spectra (low-mid-high) "green"-"gray"-"red" and "blue"-"gray"-"yellow" are used for gene.data and cpd.data respectively. The values for 'low, mid, high' can be given as color names ('red'), plot color index (2=red), and HTML-style RGB, ("\#FF0000"=red).
}
  \item{na.col}{
color used for NA's or missing values in gene.data and cpd.data. d na.col="transparent".
}
  \item{\dots}{
extra arguments passed to keggview.native or keggview.graph function.
}
}
\details{
The function KeggPathwayView is a revised version of pathview function in pathview package.
KeggPathwayView maps and renders user data on relevant pathway graphs. KeggPathwayView is a stand alone program for pathway based data integration and visualization. It also seamlessly integrates with pathway and functional analysis tools for large-scale and fully automated analysis. KeggPathwayView provides strong support for data Integration. It works with: 1) essentially all types of biological data mappable to pathways, 2) over 10 types of gene or protein IDs, and 20 types of compound or metabolite IDs, 3) pathways for over 2000 species as well as KEGG orthology, 4) varoius data attributes and formats, i.e. continuous/discrete data, matrices/vectors, single/multiple samples etc. To see mappable external gene/protein IDs do: data(gene.idtype.list), to see mappable external compound related IDs do: data(rn.list); names(rn.list). KeggPathwayView generates both native KEGG view and Graphviz views for pathways. Currently only KEGG pathways are implemented. Hopefully, pathways from Reactome, NCI and other databases will be supported in the future.
}
\value{
The result returned by KeggPathwayView function is a named list corresponding to the input pathway ids. Each element (for each pathway itself is a named list, with 2 elements ("plot.data.gene", "plot.data.cpd"). Both elements are data.frame or NULL depends on the corresponding input data gene.data and cpd.data. These data.frames record the plot data for mapped gene or compound nodes: rows are mapped genes/compounds, columns are:
\item{kegg.names }{standard KEGG IDs/Names for mapped nodes. It's Entrez Gene ID or KEGG Compound Accessions.}
\item{labels }{Node labels to be used when needed.}
\item{all.mapped }{All molecule (gene or compound) IDs mapped to this node.}
\item{type }{node type, currently 4 types are supported: "gene","enzyme", "compound" and "ortholog".}
\item{x }{x coordinate in the original KEGG pathway graph.}
\item{y }{y coordinate in the original KEGG pathway graph.}
\item{width }{node width in the original KEGG pathway graph.}
\item{height }{node height in the original KEGG pathway graph.}
\item{other columns}{columns of the mapped gene/compound data and corresponding pseudo-color codes for individual samples}
}

\author{
Wubing Zhang, Feizhen Wu, Binbin Wang
}
\seealso{
  \code{\link{pathview}}
}
\examples{

# gene_summary = BRAF_mle.gene_summary
# #Read beta score from gene summary matrix
# dd = ReadBeta(gene_summary, ctrlName = "D7_R1", treatName = "PLX7_R1")

# rownames(dd)=dd$ENTREZID
# p1 <- KeggPathwayView(gene.data  = dd[,c("Control","Treatment")],
#                   pathway.id = pathways_hsa$PathwayID[1], species="hsa")

}
