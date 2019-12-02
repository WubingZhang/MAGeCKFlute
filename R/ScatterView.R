#' Scatter plot
#'
#' Scatter plot supporting groups.
#'
#' @docType methods
#' @name ScatterView
#' @rdname ScatterView
#' @aliases ScatterView
#'
#' @param data Data frame.
#' @param x A character, specifying the x-axis.
#' @param y A character, specifying the y-axis.
#'
#' @param label An integer or a character specifying the column used as the label, default value is 0 (row names).
#' @param label.top Boolean, specifying whether label top hits.
#' @param top Integer, specifying the number of top terms in the groups to be labeled.
#' @param toplabels Character vector, specifying terms to be labeled.
#'
#' @param model One of "none" (default), "ninesquare", "volcano", and "rank".
#' @param groups Specify the colored groups. Optional groups include "top", "mid", "bottom",
#' "left", "center", "right", "topleft", "topcenter", "topright", "midleft", "midcenter",
#' "midright", "bottomleft", "bottomcenter", "bottomright".
#' @param group_col A vector of colors for specified groups.
#' @param groupnames A vector of group names to show on the legend.
#'
#' @param auto_cut Boolean, take 1.5 fold standard deviation as cutoff.
#' @param auto_cut_x Boolean, take 1.5 fold standard deviation as cutoff on x-axis.
#' @param auto_cut_y Boolean, take 1.5 fold standard deviation as cutoff on y-axis.
#' @param auto_cut_diag Boolean, take 1.5 fold standard deviation as cutoff on diagonal.
#' @param x_cut An one or two-length numeric vector, specifying the cutoff used for x-axis.
#' @param y_cut An one or two-length numeric vector, specifying the cutoff used for y-axis.
#' @param slope A numberic value indicating slope of the diagonal cutoff.
#' @param intercept A numberic value indicating intercept of the diagonal cutoff.
#'
#' @param display_cut Boolean, indicating whether display the dashed line of cutoffs.
#' @param legend Whether show the color legend.
#'
#' @param main Title of the figure.
#' @param xlab Title of x-axis
#' @param ylab Title of y-axis.
#' @param ... Other available parameters in function 'geom_text_repel'.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#' @seealso \code{\link{ScatterView}}
#'
#' @examples
#' dd = ReadBeta(mle.gene_summary)
#' ScatterView(dd, x = "dmso", y = "plx", label = "Gene",
#' x_cut = 1, y_cut = 1, groups = "topright", top = 5, display_cut = TRUE)
#'
#' @import ggplot2 ggpubr ggrepel
#' @export
#'
#'

ScatterView<-function(data, x = "x", y = "y", label = 0,
                      label.top = TRUE, top = 0, toplabels = NULL,
                      model = c("none", "ninesquare", "volcano", "rank")[1],
                      groups = NULL, group_col = NULL, groupnames = NULL,
                      auto_cut = FALSE, auto_cut_x = auto_cut,
                      auto_cut_y = auto_cut, auto_cut_diag = auto_cut,
                      x_cut = NULL, y_cut = NULL, slope = 1, intercept = NULL,
                      display_cut = FALSE, legend = FALSE,
                      # shape = 21, size = NULL, color = NULL, fill = NULL,
                      main = NULL, xlab = x, ylab = y, ...){
  requireNamespace("ggplot2", quietly=TRUE) || stop("need ggplot package")
  requireNamespace("ggrepel", quietly=TRUE) || stop("need ggrepel package")
  requireNamespace("ggpubr", quietly=TRUE) || stop("need ggpubr package")
  data = as.data.frame(data, stringsAsFactors = FALSE)
  data = na.omit(data)
  ## Add label column in the data frame.
  if(label==0) data$Label = rownames(data)
  else data$Label = as.character(data[, label])

  ## Compute the cutoff used for each dimension.
  model = tolower(model)
  if(model == "ninesquare"){
    if(length(x_cut)==0)
      x_cut = c(-CutoffCalling(data[,x], 1.5), CutoffCalling(data[,x], 1.5))
    if(length(y_cut)==0)
      y_cut = c(-CutoffCalling(data[,y], 1.5), CutoffCalling(data[,y], 1.5))
    if(length(intercept)==0)
      intercept = c(-CutoffCalling(data[,y]-data[,x], 1.5), CutoffCalling(data[,y]-data[,x], 1.5))
  }
  if(model == "volcano"){
    if(length(x_cut)==0)
      x_cut = c(-CutoffCalling(data[,x], 1.5), CutoffCalling(data[,x], 1.5))
    if(length(y_cut)==0) y_cut = -log10(0.05)
  }
  if(model == "rank"){
    if(length(x_cut)==0)
      x_cut = c(-CutoffCalling(data[,x], 1.5), CutoffCalling(data[,x], 1.5))
  }
  if(model == "none"){
    if(auto_cut_x)
      x_cut = c(-CutoffCalling(data[,x], 1.5), CutoffCalling(data[,x], 1.5))
    if(auto_cut_y)
      y_cut = c(-CutoffCalling(data[,y], 1.5), CutoffCalling(data[,y], 1.5))
    if(auto_cut_diag)
      intercept = c(-CutoffCalling(data[,y]-data[,x], 1.5), CutoffCalling(data[,y]-data[,x], 1.5))
  }
  ## Decide the colored groups
  avail_groups = c("topleft", "topright", "bottomleft", "bottomright",
                   "midleft", "topcenter", "midright", "bottomcenter", "midcenter",
                   "top", "mid", "bottom", "left", "center", "right", "none")
  ## Select the colors
  mycolour = c("#1f78b4", "#fb8072", "#33a02c", "#ff7f00",
               "#bc80bd", "#66c2a5", "#6a3d9a", "#fdb462", "#ffed6f",
               "#e78ac3", "#fdb462", "#8da0cb", "#66c2a5", "#fccde5", "#fc8d62", "#d9d9d9")
  names(mycolour) = avail_groups

  if(model == "ninesquare") groups = c("midleft", "topcenter", "midright", "bottomcenter")
  if(model == "volcano") groups = c("topleft", "topright")
  if(model == "rank") groups = c("left", "right")
  groups = intersect(groups, avail_groups)

  ## Annotate the groups in the data frame
  if(length(x_cut)>0){
    idx1 = data[,x] < min(x_cut)
    idx2 = data[,x] > max(x_cut)
  }else{
    idx1 = NA
    idx2 = NA
  }
  if(length(y_cut)>0){
    idx3 = data[,y] < min(y_cut)
    idx4 = data[,y] > max(y_cut)
  }else{
    idx3 = NA
    idx4 = NA
  }
  if(length(intercept)>0){
    idx5 = data[,y]<slope*data[,x]+min(intercept)
    idx6 = data[,y]>slope*data[,x]+max(intercept)
  }else{
    idx5 = NA; idx6 = NA
  }
  data$group="none"
  for(gr in groups){
    if(gr=="topleft") idx = cbind(idx1, idx4, idx6)
    if(gr=="topcenter") idx = cbind(!idx1, !idx2, idx4, idx6)
    if(gr=="topright") idx = cbind(idx2, idx4, idx6)
    if(gr=="midleft") idx = cbind(idx1, idx6 , !idx3, !idx4)
    if(gr=="midcenter") idx = cbind(!idx1, !idx2, !idx3, !idx4, !idx5, !idx6)
    if(gr=="midright") idx = cbind(idx2, !idx3, !idx4, idx5)
    if(gr=="bottomleft") idx = cbind(idx1, idx3, idx5)
    if(gr=="bottomcenter") idx = cbind(!idx1, !idx2, idx3, idx5)
    if(gr=="bottomright") idx = cbind(idx2, idx3, idx5)
    if(gr=="top"){
      if(length(y_cut)>0 & length(intercept)>0)
         idx = idx4 & idx6
      else if(length(y_cut)>0)
        idx = idx4
      else idx = idx6
    }
    if(gr=="mid") idx = (!idx3) & (!idx4)
    if(gr=="bottom"){
      if(length(y_cut)>0 & length(intercept)>0)
        idx = idx3 & idx5
      else if(length(y_cut)>0)
        idx = idx3
      else idx = idx5
    }
    if(gr=="left"){
      if(length(x_cut)>0 & length(intercept)>0)
        if(slope>0) idx = idx1 & idx6 else idx = idx1 & idx5
      else if(length(x_cut)>0)
        idx = idx1
      else
        if(slope>0) idx = idx6 else idx = idx5
    }
    if(gr=="center") idx = (!idx1) & (!idx2)
    if(gr=="right"){
      if(length(x_cut)>0 & length(intercept)>0)
        if(slope>0) idx = idx2 & idx5 else idx = idx2 & idx6
        else if(length(x_cut)>0)
          idx = idx2
        else
          if(slope>0) idx = idx5 else idx = idx6
    }
    ## Assign groups
    if(is.null(ncol(idx))){
      if(sum(!is.na(idx))>0) data$group[idx] = gr
      else warning("No cutpoint for group:", gr)
    }else{
      idx = idx[, !is.na(idx[1,])]
      if(is.null(ncol(idx)))
        warning("No cutpoint for group:", gr)
      else if(ncol(idx)<4 & gr=="midcenter")
        warning("No cutpoint for group:", gr)
      else
        data$group[rowSums(idx)==ncol(idx)] = gr
    }
  }
  data$group=factor(data$group, levels = unique(c(groups, "none")))
  ## Group names
  if(length(groupnames)!=length(groups)) groupnames = groups
  if(length(groups)>0) names(groupnames) = groups
  if(length(group_col)==length(groups)) mycolour[groups] = group_col
  if(length(groups)==0) mycolour["none"] = "#FF6F61"

  ## Label top gene names ##
  data$rank = top + 1
  for(g in groups){
    idx1 = data$group==g
    x_symb = 0; y_symb = 0;
    if(g=="topleft"){ x_symb = 1; y_symb = -1 }
    if(g=="topcenter"){ x_symb = 0; y_symb = -1 }
    if(g=="topright"){ x_symb = -1; y_symb = -1 }
    if(g=="midleft"){ x_symb = 1; y_symb = 0 }
    if(g=="midright"){ x_symb = -1; y_symb = 0 }
    if(g=="bottomleft"){ x_symb = 1; y_symb = 1 }
    if(g=="bottomcenter"){ x_symb = 0; y_symb = 1 }
    if(g=="bottomright"){ x_symb = -1; y_symb = 1 }
    if(g=="top"){ x_symb = 0; y_symb = -1 }
    if(g=="bottom"){ x_symb = 0; y_symb = 1 }
    if(g=="left"){ x_symb = 1; y_symb = 0 }
    if(g=="right"){ x_symb = -1; y_symb = 0 }
    tmp = data[,c(x,y)]
    tmp[,x] = (tmp[,x]-min(tmp[,x])) / (max(tmp[,x])-min(tmp[,x]))
    tmp[,y] = (tmp[,y]-min(tmp[,y])) / (max(tmp[,y])-min(tmp[,y]))
    data$rank[idx1] = rank((x_symb*tmp[,x]+y_symb*tmp[,y])[idx1])
  }
  data$rank[data$rank==0] = Inf
  if(mode(toplabels)=="list"){
    data$Label[data$rank>top & !(data$Label %in% unlist(toplabels))] = ""
    data$Color = data$Label;
    if(length(toplabels)>0){
      tmp = stack(toplabels)
      tmp = tmp[!duplicated(tmp[,1]), ]
      rownames(tmp) = tmp[,1]
      data$Color[data$Color%in%tmp[,1]] = as.character(tmp[data$Color[data$Color%in%tmp[,1]], 2])
      data$Color[!(data$Color%in%tmp[,2]) & data$Color!=""] = "Top hits"
    }
  }else{
    data$Label[data$rank>top & !(data$Label %in% toplabels)] = ""
  }
  ## Plot the scatter figure ##
  gg = data
  if(mode(toplabels)!="list"){
    p = ggplot(gg, aes_string(x, y, label="Label", color = "group"), alpha = 0.8)
    p = p + geom_point(shape = 16)
    p = p + scale_color_manual(values = mycolour, labels = groupnames)
  }else{
    p = ggplot(gg, aes_string(x, y, label="Label"), alpha = 0.8)
    p = p + geom_point(shape = 16, color = "#d9d9d9")
    p = p + geom_point(aes(color = Color), shape = 16, data = gg[gg$Color!="", ])
  }
  if(label.top)
    p = p + ggrepel::geom_text_repel(...)
  if(display_cut){
    if(length(x_cut)>0)
      p = p + geom_vline(xintercept = x_cut,linetype = "dotted")
    if(length(y_cut)>0)
      p = p + geom_hline(yintercept = y_cut,linetype = "dotted")
    if(length(intercept)>0)
      p = p + geom_abline(slope=slope, intercept=intercept, linetype = "dotted")
  }
  p = p + labs(x=xlab, y = ylab, title = main, color = NULL)
  p = p + ggpubr::theme_pubr() + theme(plot.title = element_text(hjust = 0.5))
  if(!legend) p = p + theme(legend.position = "none")
  return(p)
}


#' Quantile of normal distribution.
#'
#' Compute cutoff from a normal-distributed vector.
#'
#' @docType methods
#' @name CutoffCalling
#' @rdname CutoffCalling
#'
#' @param d A numeric vector.
#' @param scale Boolean or numeric, specifying how many standard deviation will be used as cutoff.
#'
#' @return A numeric value.
#' @export
#' @examples
#' CutoffCalling(rnorm(10000))

CutoffCalling=function(d, scale=1){
  param=1
  if(is.logical(scale) & scale){
    param = round(length(d) / 20000, digits = 1)
  }else if(is.numeric(scale)){param = scale}

  Control_mean=0
  sorted_beta=sort(abs(d))
  temp=quantile(sorted_beta,0.68)
  temp_2=qnorm(0.84)
  cutoff=round(temp/temp_2,digits = 3)
  names(cutoff)=NULL
  cutoff=cutoff*param
  return(cutoff)
}
