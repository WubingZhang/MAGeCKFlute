#' Call cutoff
#'
#' Calculate standard deviation as cutoff for a numeric vector
#'
#' @docType methods
#' @name CutoffCalling
#' @rdname CutoffCalling
#'
#' @param d A numeric vector.
#' @param scale Boolean or numeric, whether scale cutoff to whole genome level,
#' or how many standard deviation will be used as cutoff.
#'
#' @return A numeric value.
#' @examples
#' CutoffCalling(rnorm(10000))
#' @export

CutoffCalling=function(d, scale=FALSE){
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
