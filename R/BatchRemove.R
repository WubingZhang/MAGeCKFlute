#' Batch effect removal
#'
#' @docType methods
#' @name BatchRemove
#' @rdname BatchRemove
#'
#' @param mat A data frame, each row is a gene, and each column is a sample.
#' @param batchMat A data frame, the first column should be `Samples`(matched colnames of mat)
#' and the second column is `Batch`. The remaining columns could be Covariates.
#' @param log2trans Boolean, specifying whether do logarithmic transformation before batch removal.
#'
#' @return A list contrains two objects, including \code{data} and \code{p}.
#'
#' @author Wubing Zhang
#'
#' @seealso \code{\link[sva]{ComBat}}
#'
#' @examples
#' edata = matrix(c(rnorm(2000, 5), rnorm(2000, 8)), 1000)
#' colnames(edata) = paste0("s", 1:4)
#' batchMat = data.frame(sample = colnames(edata), batch = rep(1:2, each = 2))
#' edata1 = BatchRemove(edata, batchMat)
#' print(edata1$p)
#'
#' @importFrom sva ComBat
#'
#' @export
#'
BatchRemove <- function(mat, batchMat, log2trans=FALSE){
  requireNamespace("sva", quietly=TRUE) || stop("need sva package")
  mat = as.data.frame(mat, stringsAsFactors = FALSE)
  batchMat = as.data.frame(batchMat, stringsAsFactors = FALSE)
  colnames(batchMat)[1:2] = c("Sample", "Batch")
  rownames(batchMat) = batchMat[,1]
  batch = batchMat
  requireNamespace("sva")
  # mat = as.matrix(mat)
  # if(mode(mat)!="numeric") stop("Numeric data matrix is needed!")
  # batch = as.matrix(batchMat)
  # rownames(batch) = batch[,1]

  ## load batch matrix
  index=intersect(batchMat[,1], colnames(mat))
  if(length(index)<2) stop("Too less samples ...")
  dt = mat[, index]
  if(log2trans) dt = log(dt+1)
  ##
  tmp=dt
  var0=numeric()
  for (i in unique(batch[,2])){
    temp=as.character(batch[batch[,2]==i,1])
    if (length(temp)>=2){
      temp2=tmp[,temp]
      var.data.tmp <- apply(temp2, 1,function(x) var(x))
      var.data.tmp0=which(var.data.tmp==0)
      var0=c(var0,var.data.tmp0)
    }
  }
  var0 = unique(var0)
  tmp2 = tmp
  if(length(var0)>0) tmp2=tmp[setdiff(1:nrow(tmp),var0),]
  if(ncol(batch)>2){
    mod = as.data.frame(batch[,3:ncol(batch)])
    mod = model.matrix(~., data=mod)
  }else mod = NULL
  if(length(unique(batch[,2]))<2){
    res = tmp2
  }else{res <- sva::ComBat(as.matrix(tmp2), batch = batch[,2], mod = mod)}
  # if(positive) res[res<0] = 0

  if(length(setdiff(colnames(mat), colnames(res)))>0){
    clnames = setdiff(colnames(mat),colnames(res))
    res=cbind.data.frame(mat[setdiff(1:nrow(tmp),var0), clnames], as.data.frame(res))
    colnames(res)[1:length(clnames)] = clnames
  }
  if(length(var0)>0){
    res = rbind.data.frame(as.data.frame(res), mat[var0, colnames(res)])
  }

  # dt2 = matrix(as.numeric(res[,index]), ncol = length(index))
  # colnames(dt2) = index
  # rownames(dt2) = rownames(res)
  dt2 = res[,index]
  pca1 = NULL
  pca2 = NULL
  p1 = NULL

  pca1 = prcomp(t(dt))$x[,1:2]
  pca2 = prcomp(t(dt2))$x[,1:2]

  gg1 = as.data.frame(pca1, stringsAsFactors=FALSE)
  gg1$col = as.character(batch[rownames(gg1),2])
  gg1$group = factor("Before batch removal", "Before batch removal")
  gg2 = as.data.frame(pca2,stringsAsFactors=FALSE)
  gg2$col = as.character(batch[rownames(gg2),2])
  gg2$group = factor("After batch removal", "After batch removal")
  if(ncol(batch)>2){
    gg1$shape = as.character(batch[rownames(gg1),3])
    gg2$shape = batch[rownames(gg2),3]
  }else{ gg1$shape = "NA"; gg2$shape = "NA" }
  gg = rbind.data.frame(gg1, gg2)

  #====plot PC1 and PC2=====
  p1 = ggplot(gg)
  p1 = p1 + geom_point(aes(x=PC1, y=PC2, color=col, shape=shape),size = 1)
  p1 = p1 + scale_color_discrete(name="Batch", breaks = unique(gg$col))
  p1 = p1 + scale_shape_discrete(name="Batch", breaks = unique(gg$col))
  p1 = p1 + theme_bw(14)+theme(plot.title = element_text(hjust = 0.5,size=12))
  p1 = p1 + facet_grid(~group, switch = "y", scales="free")
  p1 = p1 + theme(legend.title=element_blank())
  # ggsave(file.path(outdir, paste0(prefix, "_PCA_BatchRemoval.png")), p1, width = 10, height = 4)

  ##====Clustering========
  # if(cluster){
  #   filename = file.path(outdir, paste0(prefix, "_hclust_batchremoval.pdf"))
  #   if(is.na(hclust.height) | is.na(hclust.width))
  #     pdf(filename, height=0.5*nrow(batch)+2, width=0.05*nrow(batch)+3)
  #   else
  #     pdf(filename, height=hclust.height, width=hclust.width)
  #   cc2 = cor(dt)
  #   cc=cor(dt2)
  #   if(ncol(batch)>2){
  #     hclustView(cc2, label_cols = batch[rownames(cc2),3], bar_cols = batch[rownames(cc2),2:ncol(batch)], main = "Before batch removal")
  #     hclustView(cc, label_cols = batch[rownames(cc),3], bar_cols = batch[rownames(cc),2:ncol(batch)], main = "After batch removal")
  #
  #   }else{
  #     hclustView(cc2, label_cols = batch[rownames(cc2),2], main = "Before batch removal")
  #     hclustView(cc, label_cols = batch[rownames(cc),2], main = "After batch removal")
  #   }
  #   dev.off()
  # }
  return(list(data=res, p=p1))
}
