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
#' @param batchMat data frame or matrix, which has three columns of Samples matched colname of mat, Batch, and Covariate.
#' @param cluster boolean, specifing whether do cluster analysis before and after batch removal
#' @param log2trans boolean, specifing whether do log2 transition before batch removal
#' @param ... other parameters in pheatmap
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
#' res = BatchRemove(beta, batchMat, cluster = FALSE)
#' @importFrom sva ComBat
#'
#' @export
#'

BatchRemove <- function(mat, batchMat, cluster=TRUE, log2trans=FALSE, ...){
  requireNamespace("sva")
  batch = batchMat
  ## load batch matrix
  index=match(batch[,1],colnames(mat))
  dt=as.matrix(mat[,index])
  ##
  tmp=dt
  var0=numeric()
  for (i in unique(batch[,2])){
    temp=as.character(batch[batch[,2]==i,1])
    if (length(temp)>=2){
      temp2=tmp[,temp]
      var.data.tmp <- apply(temp2, 1,function(x) var(x))
      var.data.tmp0=var.data.tmp[var.data.tmp==0]
      var0=c(var0,var.data.tmp0)
    }
  }
  if(length(var0)>0) tmp2=tmp[setdiff(rownames(tmp),names(var0)),]
  if(log2trans) tmp2=log2(tmp2+0.000001)    #====Revised
  res <- sva::ComBat(tmp2, batch = batch[,2], mod = batch[,3])
  if(log2trans) res = 2^res  ##revised

  if(length(setdiff(colnames(mat), colnames(res)))>0)
    res2=cbind(mat[rownames(res), setdiff(colnames(mat),colnames(res))], res)
  if(length(setdiff(rownames(mat), rownames(res2)))>0)
    after = rbind(res2, mat[setdiff(rownames(mat),rownames(res2)),colnames(res2)])

  ##====Clustering========
  if(cluster){
    display_numbers = TRUE
    if(ncol(dt)>10){display_numbers = FALSE}
    HeatmapView(res2[, colnames(dt)], filename = 'ClusterAfterBatchRemove.pdf', display_numbers = display_numbers, ...)
    HeatmapView(dt, filename = 'ClusterBeforeBatchRemove.pdf', display_numbers = display_numbers, ...)
  }
  return(after_removal=after)
}
