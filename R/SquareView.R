#' Scatter plot of 9-Square
#'
#' Plot a scatter plot with Control beta score as x-axis and Treatment beta score as y-axis,
#' and colored treatment related genes.
#'
#' @docType methods
#' @name SquareView
#' @rdname SquareView
#' @aliases squareview
#'
#' @param beta data frame, which has columns of 'Gene', \code{ctrlname} and \code{treatname}.
#' @param ctrlname a character, specifying the names of control samples.
#' @param treatname a character, specifying the name of treatment samples.
#' @param genelist character vector, specifying labeled genes.
#' @param scale_cutoff boolean or numeric, whether scale cutoff to whole genome level,
#' or how many standard deviation will be used as cutoff.
#' @param main as in 'plot'.
#' @param filename figure file name to create on disk. Default filename="NULL", which means
#' don't save the figure on disk.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#' @note See the vignette for an example of SquareView.
#' Note that the source code of \code{SquareView} is very simple.
#' The source can be found by typing \code{MAGeCKFlute:::SquareView}
#' or \code{getMethod("SquareView")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/SquareView.R}
#' Users should find it easy to customize this function.
#'
#' @seealso \code{\link{DensityDiffView}}   \code{\link{DensityView}}
#' @seealso \code{\link{ViolinView}}   \code{\link{MAView}}
#' @seealso \code{\link{CellCycleView}}  \code{\link{EnrichedView}}
#' @seealso \code{\link{EnrichedGSEView}}  \code{\link{KeggPathwayView}}
#' @seealso \code{\link{RankView}}    \code{\link{ScatterView}}
#'
#' @examples
#' data(MLE_Data)
#' # Read beta score from gene summary table in MAGeCK MLE results
#' dd = ReadBeta(MLE_Data, organism="hsa")
#' SquareView(dd, ctrlname = "D7_R1", treatname = "PLX7_R1")
#'
#'
#' @importFrom ggExtra ggMarginal
#' @importFrom ggrepel geom_text_repel
#' @export
#'
#'

#Plot square
SquareView<-function(beta, ctrlname="Control",treatname="Treatment", genelist=c(),
                     scale_cutoff=1.2, main=NULL, filename=NULL){
  requireNamespace("ggExtra", quietly=TRUE) || stop("need ggExtra package")
  requireNamespace("ggrepel", quietly=TRUE) || stop("need ggrepel package")

  loginfo(paste("Square plot for",main, "beta scores ..."))
  beta$group="Others"
  beta$text = beta$Gene
  beta$Control = rowMeans(beta[, ctrlname, drop= FALSE])
  beta$Treatment = rowMeans(beta[, treatname, drop= FALSE])
  Control_cutoff = CutoffCalling(beta$Control, scale=scale_cutoff)
  drug_cutoff = CutoffCalling(beta$Treatment, scale=scale_cutoff)
  x=beta$Control; y=beta$Treatment
  intercept=CutoffCalling(y-x,scale=scale_cutoff)
  y_max=x+intercept; y_min=x-intercept
  idx0=y>y_max | y<y_min

  idx1 = -Control_cutoff<beta$Control
  idx2 = beta$Control<Control_cutoff
  idx3 = beta$Treatment>drug_cutoff
  idx=idx0&idx1&idx2&idx3
  beta$group[idx]="Group2" #proliferation
  beta[idx,] = beta[which(idx)[order(beta$Treatment[idx], decreasing = TRUE)],]
  beta$text[which(idx)[rank(beta$Treatment[idx])<(length(beta$Treatment[idx])-9)]] = NA
  beta[idx,] = beta[which(idx)[order(beta$Treatment[idx], decreasing = TRUE)],]
  idx3=beta$Treatment<(-drug_cutoff)
  idx=idx0&idx1&idx2&idx3
  beta$group[idx]="Group4" #anti-proliferation
  beta$text[which(idx)[rank(beta$Treatment[idx])>10]] = NA
  beta[idx,] = beta[which(idx)[order(beta$Treatment[idx], decreasing = FALSE)],]
  idx1=-drug_cutoff<beta$Treatment
  idx2=beta$Treatment<drug_cutoff
  idx3=beta$Control>Control_cutoff
  idx=idx0&idx1&idx2&idx3
  beta$group[idx]="Group3" #proliferation
  beta$text[which(idx)[rank(beta$Control[idx])<(length(beta$Control[idx])-9)]] = NA
  beta[idx,] = beta[which(idx)[order(beta$Control[idx], decreasing = TRUE)],]
  idx3=beta$Control<(-Control_cutoff)
  idx=idx0&idx1&idx2&idx3
  beta$group[idx]="Group1" #anti-proliferation
  beta$text[which(idx)[rank(beta$Control[idx])>10]] = NA
  beta[idx,] = beta[which(idx)[order(beta$Control[idx], decreasing = FALSE)],]
  #===================
  x_min=round(min(beta$Control),2)+0.02
  x_max=round(max(beta$Control),2)-0.02
  y_min=round(min(beta$Treatment),2)+0.02
  y_max=round(max(beta$Treatment),2)-0.02

  gg = beta
  gg$group=factor(gg$group,levels = c("Group1","Group2","Group3","Group4","Others"))
  #===============
  # gg=gg[,c("Gene","Treatment","Control","group","ENTREZID")]
  mycolour=c("Others"="aliceblue",  "Group2"="#ff7f00", "Group3"="#ffff33",
             "Group4"="#984ea3", "Group1"="#4daf4a" )
  idx1 = gg$Gene %in% genelist
  idx2 = !(gg$group=="Others" | is.na(gg$text))
  label_gg = gg[idx1|idx2,]
  col_label = rep("DeepSkyBlue4",nrow(label_gg))
  col_label[label_gg$group=="Others"]="gray60"
  p=ggplot(gg,aes(x=Control,y=Treatment,colour=group,fill=group))
  p=p+geom_point(shape=".",alpha=1/1,size = 1)+scale_color_manual(values=mycolour)
  p=p+geom_jitter(size = 1)
  p=p+theme_bw(14)+theme(plot.title = element_text(hjust = 0.5,size=12))
  p=p+geom_vline(xintercept = Control_cutoff,linetype = "dotted")
  p=p+geom_vline(xintercept = -Control_cutoff,linetype = "dotted")
  p=p+geom_hline(yintercept = drug_cutoff,linetype = "dotted")
  p=p+geom_hline(yintercept = -drug_cutoff,linetype = "dotted")
  p=p+geom_abline(slope=1, intercept=+intercept,linetype = "dotted")
  p=p+geom_abline(slope=1, intercept=-intercept,linetype = "dotted")
  p=p+labs(x=paste0(ctrlname, collapse = "/"), y=paste0(treatname, collapse = "/"),title=main)
  p=p+theme(legend.position="bottom")+theme(legend.title=element_blank())
  p=p+guides(col = guide_legend(ncol = 3,byrow = TRUE))
  p=p+ggrepel::geom_text_repel(aes(x=Control,y=Treatment,label=Gene), color=col_label, data=label_gg)
  # p=p+xlim(x_min,x_max)+ylim(y_min,y_max)
  p=p+annotate("text",color="red",
               label=paste("",as.character(dim(beta[beta$group=="Group2",])[1]),sep=""),
               x=0, y=y_max,vjust=0,hjust = 0.5)
  p=p+annotate("text",color="red",
               label=paste("",as.character(dim(beta[beta$group=="Group4",])[1]),sep=""),
               x=0, y=y_min,vjust=1,hjust = 0.5)
  p=p+annotate("text",color="red",
               label=paste("",as.character(dim(beta[beta$group=="Group3",])[1]),sep=""),
               x=x_max, y=0,vjust=0.5,hjust = 1)
  p=p+annotate("text",color="red",
               label=paste("",as.character(dim(beta[beta$group=="Group1",])[1]),sep=""),
               x=x_min, y=0,vjust=0.5,hjust = 0)
  p=suppressWarnings(ggExtra::ggMarginal(p, type="histogram",bins=50))
  p$data = gg[, c("Gene", "ENTREZID", ctrlname, treatname, "group")]
  p$data = p$data[order(p$data$group),]
  # grid.arrange(p,ncol = 1)
  if(!is.null(filename)){
      write.table(p$data, gsub("\\....$", ".txt", filename),
                sep = "\t", quote = FALSE, row.names = FALSE)
      ggsave(filename, p, units="in",width=10,height=8)
  }
  return(p)
}
