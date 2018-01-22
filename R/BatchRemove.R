#' Batch effect removal
#'
#' Remove batch effect
#'
#' @docType methods
#' @name BatchRemove
#' @rdname BatchRemove
#' @aliases batchremove
#'
#' @param mat data frame or matrix, in which each column represents one sample.
#' @param batchMat data frame or matrix, which has at least three columns, including Samples matched colname of mat, Batch, and Covariates.
#' @param cov boolean, specifying if provide mod in ComBat
#' @param log2trans boolean, specifying whether do log2 transition before batch removal
#' @param pca boolean, specifying whether do principle component analysis before and after batch removal
#' @param cluster boolean, specifing whether do cluster analysis before and after batch removal
#' @param prefix character, specifying prefix of output figures, only needed if cluster/pca is TRUE.
#' @param outdir  output directory on disk
#'
#' @return matrix of data after batch removal.
#'
#' @author Wubing Zhang
#'
#' @note See the vignette for an example of BatchRemove
#' Note that the source code of \code{BatchRemove} is very simple.
#' The source can be found by typing \code{MAGeCKFlute:::BatchRemove}
#' or \code{getMethod("BatchRemove")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/BatchRemove.R}
#' Users should find it easy to customize this function.
#'
#' @seealso \code{\link[sva]{ComBat}}
#'
#' @examples
#' beta = ReadBeta(MLE_Data, organism="hsa")
#' batchMat = data.frame(samples = c("D7_R1", "D7_R2", "PLX7_R1", "PLX7_R2"),
#'                       batch = c(1,2,1,2), cov = c(1,1,2,2))
#' res = BatchRemove(beta, batchMat, prefix="Test")
#'
#' @importFrom sva ComBat
#' @importFrom sva ComBat
#'
#' @export
#'

BatchRemove <- function(mat, batchMat, cov=FALSE, log2trans=FALSE, pca=TRUE, cluster=TRUE, prefix=NA, outdir="."){
  requireNamespace("sva")
  batch = as.matrix(batchMat)
  rownames(batch) = batch[,1]
  ## load batch matrix
  index=intersect(batch[,1], colnames(mat))
  dt=as.matrix(mat[,index])
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
  tmp2 = tmp
  if(length(var0)>0) tmp2=tmp[setdiff(1:nrow(tmp2),var0),]
  if(cov){
    mod = as.data.frame(batch[,3:ncol(batch)])
    mod = model.matrix(~., data=mod)
  }else mod = NULL
  if(length(unique(batch[,2]))<2){
    res = tmp2
  }else{res <- sva::ComBat(tmp2, batch = batch[,2], mod = mod)}
  if(length(setdiff(colnames(mat), colnames(res)))>0)
    res=cbind(mat[rownames(res), setdiff(colnames(mat),colnames(res))], res)
  if(length(setdiff(rownames(mat), rownames(res)))>0)
    res = rbind(res, mat[setdiff(rownames(mat),rownames(res)),colnames(res)])

  dt2 = as.matrix(res[,index])
  pca1 = NULL
  pca2 = NULL
  p1 = NULL
  if(pca){
    pca1 = prcomp(t(dt))$x[,1:2]
    pca2 = prcomp(t(dt2))$x[,1:2]

    gg1 = as.data.frame(pca1, stringsAsFactors=FALSE)
    gg1$col = batch[rownames(gg1),2]
    gg1$group = factor("Before batch removal", "Before batch removal")
    gg2 = as.data.frame(pca2,stringsAsFactors=FALSE)
    gg2$col = batch[rownames(gg2),2]
    gg2$group = factor("After batch removal", "After batch removal")
    if(ncol(batch)>2){
      gg1$shape = as.character(batch[rownames(gg1),3])
      gg2$shape = batch[rownames(gg2),3]
    }else{ gg1$shape = "NA"; gg2$shape = "NA" }
    gg = rbind.data.frame(gg1, gg2)

    #====plot PC1 and PC2=====
    p1 = ggplot(gg)
    p1 = p1 + geom_point(aes(x=PC1, y=PC2, color=col, shape=shape),size = 1)
    p1 = p1 + theme_bw(14)+theme(plot.title = element_text(hjust = 0.5,size=12))
    p1 = p1 + facet_grid(~group, switch = "y", scales="free")
    p1 = p1 + theme(legend.title=element_blank())
    ggsave(file.path(outdir, paste0(prefix, "_PCA_BatchRemoval.png")), p1, width = 10, height = 4)
  }

  ##====Clustering========
  if(cluster){
    filename = file.path(outdir, paste0(prefix, "_hclust_batchremoval.pdf"))

    pdf(filename, height=0.2*nrow(batch)+2, width=0.02*nrow(batch)+3)
    cc2 = cor(dt)
    cc=cor(dt2)
    if(ncol(batch)>2){
      hclustView(cc2, label_cols = batch[rownames(cc2),3], bar_cols = batch[rownames(cc2),2:ncol(batch)], main = "Before batch removal")
      hclustView(cc, label_cols = batch[rownames(cc),3], bar_cols = batch[rownames(cc),2:ncol(batch)], main = "After batch removal")

    }else{
      hclustView(cc2, label_cols = batch[rownames(cc2),2], main = "Before batch removal")
      hclustView(cc, label_cols = batch[rownames(cc),2], main = "After batch removal")
    }
    dev.off()

    loginfo("DONE!")
  }
  return(list(data=res, pca1=pca1, pca2=pca2, p=p1))
}
