#' Call cutoff
#'
#' Calculate standard deviation as cutoff for a numeric vector
#'
#' @docType methods
#' @name CutoffCalling
#' @rdname CutoffCalling
#'
#' @param d A numeric vector.
#' @param scale Boolean or numeric, whether scale cutoff to whole genome level, or how many standard deviation will be used as cutoff.
#'
#' @return A numeric value.
#'
#' @author Wubing Zhang
#'
#' @note See the vignette for an example of cutoff calling.
#' Note that the source code of \code{CutoffCalling} is very simple.
#' The source can be found by typing \code{MAGeCKFlute:::CutoffCalling}
#' or \code{getMethod("CutoffCalling")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/CutoffCalling.R}
#' Users should find it easy to customize this function.
#'
#'

CutoffCalling=function(d, scale=FALSE){
  param=1
  if(class(scale)=="logical" & scale){
    param=round(length(d) / 18000,digits = 1)
  }else if(class(scale)=="numeric"){param = scale}

  Control_mean=0
  sorted_beta=sort(abs(d))
  temp=quantile(sorted_beta,0.68)
  temp_2=qnorm(0.84)
  cutoff=round(temp/temp_2,digits = 3)
  names(cutoff)=NULL
  cutoff=cutoff*param
  if(cutoff==0){
    stop("Cutoff can not be zero!")
  }
  return(cutoff)
}
