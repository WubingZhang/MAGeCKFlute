#' Map values to colors
#'
#' Map values to colors
#'
#' @docType methods
#' @name getCols
#' @rdname getCols
#'
#' @param x a numeric vector
#' @param palette an integer
#'
#' @return a vector of colors corresponding to input vector.
#'
#' @author Wubing Zhang
#'
#'
#' @importFrom colorspace rainbow_hcl
#' @importFrom scales gradient_n_pal
#' @importFrom scales brewer_pal
#'
getCols <- function(x, palette=1){
  requireNamespace("scales")
  if(palette==1) cols <- colorspace::rainbow_hcl(length(unique(x)), c = 70, l  = 50)
  if(palette>1) cols = scales::gradient_n_pal(scales::brewer_pal("div", palette-1)(8))(seq(0, 1, length.out = length(unique(x))))

  names(cols) = unique(x)
  return(cols[x])
}
