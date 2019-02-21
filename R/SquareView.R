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
#' @param beta Data frame, including columns of \code{ctrlname} and \code{treatname}, with Gene Symbol as rowname.
#' @param ctrlname A character, specifying the names of control samples.
#' @param treatname A character, specifying the name of treatment samples.
#' @param label An integer or a character specifying the column used as the label, default value is 0 (row names).
#' @param label.top Boolean, whether label the top selected genes, default label the top 10 genes in each group.
#' @param top Integer, specifying the number of top selected genes to be labeled. Default is 5.
#' @param genelist Character vector, specifying labeled genes.
#' @param scale_cutoff Numeric, specifying the number of standard deviation to be used as cutoff.
#' @param main As in 'plot'.
#' @param filename Figure file name to create on disk. Default filename="NULL", which means
#' don't save the figure on disk.
#' @param width As in ggsave.
#' @param height As in ggsave.
#' @param ... Other available parameters in function 'ggsave'.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#' @seealso \code{\link{ScatterView}}
#'
#' @examples
#' data(mle.gene_summary)
#' # Read beta score from gene summary table in MAGeCK MLE results
#' dd = ReadBeta(mle.gene_summary, organism="hsa")
#' SquareView(dd, ctrlname = "dmso", treatname = "plx", label = "Gene")
#'
#'
#' @importFrom ggExtra ggMarginal
#' @importFrom ggrepel geom_text_repel
#' @export
#'
#'

#Plot square
SquareView<-function(beta, ctrlname = "Control", treatname = "Treatment",
                     label = 0, label.top = TRUE, top=5, genelist=c(),
                     scale_cutoff = 1.96, main=NULL,
                     filename=NULL, width=5, height=4, ...){
  requireNamespace("ggExtra", quietly=TRUE) || stop("need ggExtra package")
  requireNamespace("ggrepel", quietly=TRUE) || stop("need ggrepel package")

  message(Sys.time(), " # Square plot for ", main, " beta scores ...")
  if(!all(c(ctrlname, treatname) %in% colnames(beta))){
    stop("No sample found in the data.")
  }
  if(label==0)
    beta$Gene = rownames(beta)
  else
    beta$Gene = beta[, label]
  beta$group="Others"
  beta$text = beta$Gene
  beta$Control = rowMeans(beta[, ctrlname, drop= FALSE])
  beta$Treatment = rowMeans(beta[, treatname, drop= FALSE])
  Control_cutoff = CutoffCalling(beta$Control, scale=scale_cutoff)
  drug_cutoff = CutoffCalling(beta$Treatment, scale=scale_cutoff)
  x=beta$Control; y=beta$Treatment
  beta$diff = y-x
  intercept=CutoffCalling(y-x, scale = 1)
  y_max=x+intercept; y_min=x-intercept
  idx0=y>y_max | y<y_min

  beta$rank = nrow(beta)
  idx1 = -Control_cutoff<beta$Control
  idx2 = beta$Control<Control_cutoff
  idx3 = beta$Treatment>drug_cutoff
  idx=idx0&idx1&idx2&idx3
  beta$group[idx]="Group2" #proliferation
  beta$rank[idx] = sum(idx) - rank(beta$Treatment[idx]) + 1

  idx3=beta$Treatment<(-drug_cutoff)
  idx=idx0&idx1&idx2&idx3
  beta$group[idx]="Group4" #anti-proliferation
  beta$rank[idx] = rank(beta$Treatment[idx])

  idx1=-drug_cutoff<beta$Treatment
  idx2=beta$Treatment<drug_cutoff
  idx3=beta$Control>Control_cutoff
  idx=idx0&idx1&idx2&idx3
  beta$group[idx]="Group3" #proliferation
  beta$rank[idx] = sum(idx) - rank(beta$Control[idx]) + 1

  idx3=beta$Control<(-Control_cutoff)
  idx=idx0&idx1&idx2&idx3
  beta$group[idx]="Group1" #anti-proliferation
  beta$rank[idx] = rank(beta$Control[idx])
  #===================
  beta$text[beta$group=="Others"] = NA
  beta$text[beta$rank > top] = NA
  tmp = ifelse(label.top, 0.1, 0.4)
  x_min = round(min(beta$Control[beta$group != "Others"]),2) - tmp
  x_max = round(max(beta$Control[beta$group != "Others"]),2) + tmp
  y_min = round(min(beta$Treatment[beta$group != "Others"]),2) - tmp
  y_max = round(max(beta$Treatment[beta$group != "Others"]),2) + tmp

  idx1 = (x_min<=beta$Control & beta$Control<=x_max)
  idx2 = (y_min<=beta$Treatment  & beta$Treatment<=y_max)
  idx=idx1&idx2
  gg = beta[idx,]
  gg$group=factor(gg$group, levels = c("Group1","Group2","Group3","Group4","Others"))
  #===============
  mycolour=c("Others"="gray80",  "Group2"="#ff7f00", "Group3"="#005CB7",
             "Group4"="#984ea3", "Group1"="#4daf4a" )
  idx1 = gg$Gene %in% genelist
  idx2 = !(gg$group=="Others" | is.na(gg$text))
  label_gg = gg[idx1|idx2,]
  col_label = rep("#004b84",nrow(label_gg))
  col_label[label_gg$group=="Others"]="gray60"
  p = ggplot(gg,aes(x=Control,y=Treatment,colour=group,fill=group))
  p = p + geom_point(shape=".", alpha=0.8, size = 1)+scale_color_manual(values=mycolour)
  p = p + geom_jitter(size = 1)
  p = p + geom_vline(xintercept = Control_cutoff,linetype = "dotted")
  p = p + geom_vline(xintercept = -Control_cutoff,linetype = "dotted")
  p = p + geom_hline(yintercept = drug_cutoff,linetype = "dotted")
  p = p + geom_hline(yintercept = -drug_cutoff,linetype = "dotted")
  p = p + geom_abline(slope=1, intercept=+intercept,linetype = "dotted")
  p = p + geom_abline(slope=1, intercept=-intercept,linetype = "dotted")
  p = p + labs(x="Control beta score", y = "Treatment beta score")
  # p=p+labs(x=paste0(ctrlname, collapse = " & "), y=paste0(treatname, collapse = " & "), title=main)
  p=p+guides(col = guide_legend(ncol = 3, byrow = TRUE))
  if(label.top)
    p = p + ggrepel::geom_text_repel(aes(x=Control,y=Treatment,label=Gene), color=col_label, data=label_gg)
  p=p+annotate("text",color="red",
               label=paste("",as.character(dim(beta[beta$group=="Group2",])[1]),sep=""),
               x = -Control_cutoff, y=y_max,vjust=0, hjust = 0)
  p=p+annotate("text", color="red",
               label=paste("",as.character(dim(beta[beta$group=="Group4",])[1]),sep=""),
               x = Control_cutoff, y=y_min,vjust=1,hjust = 1)
  p=p+annotate("text",color="red",
               label=paste("",as.character(dim(beta[beta$group=="Group3",])[1]),sep=""),
               x = x_max, y = -drug_cutoff, vjust=0, hjust = 1)
  p=p+annotate("text",color="red",
               label=paste("",as.character(dim(beta[beta$group=="Group1",])[1]),sep=""),
               x=x_min, y=drug_cutoff,vjust=1,hjust = 0)
  p = p + theme(legend.key = element_rect(fill = "transparent"))
  p = p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank())
  p=p+theme(legend.position="none")+theme(legend.title=element_blank())

  # p = suppressWarnings(ggExtra::ggMarginal(p, type="histogram", bins=50, fill = "gray80"))
  p$data = beta
  p$data = p$data[order(p$data$group), ]

  if(!is.null(filename)){
      write.table(p$data, gsub("\\....$", ".txt", filename),
                sep = "\t", quote = FALSE, row.names = FALSE)
      ggsave(filename, p, units="in", width=width, height=height, ...)
  }
  return(p)
}
