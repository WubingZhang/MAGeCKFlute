#' Normalize gene beta scores
#'
#' Two normalization methods are available. \code{cell_cycle} method normalizes gene beta scores
#'  based on positive control genes in CRISPR screening. \code{loess} method normalizes gene
#'  beta scores using loess.
#'
#' @docType methods
#' @name NormalizeBeta
#' @rdname NormalizeBeta
#' @aliases normalizebeta
#'
#' @param beta Data frame.
#' @param id An integer specifying the column of gene.
#' @param method Character, one of 'cell_cycle'(default) and 'loess'.
#' or character string giving the name of the table column containing the gene names.
#' @param posControl A character vector, specifying a list of positive control genes.
#' @param samples Character vector, specifying the sample names in \emph{beta} columns.
#' If NULL (default), take all \emph{beta} columns as samples.
#'
#' @return A data frame with same format as input data \emph{beta}.
#'
#' @details In CRISPR screens, cells treated with different conditions (e.g., with or without drug)
#' may have different proliferation rates. So it's necessary to normalize the proliferation rate
#' based on defined positive control genes among samples. After normalization, the beta scores are
#' comparable across samples. \code{loess} is another optional normalization method, which is used
#' to normalize array data before.
#'
#' @author Wubing Zhang
#'
#' @examples
#' data(mle.gene_summary)
#' data(Zuber_Essential)
#' # Read beta score from gene summary table in MAGeCK MLE results
#' dd = ReadBeta(mle.gene_summary)
#' #Cell Cycle normalization
#' dd_essential = NormalizeBeta(dd, samples=c("dmso", "plx"),
#'     method="cell_cycle", posControl = Zuber_Essential$GeneSymbol)
#' head(dd_essential)
#'
#' #Optional loess normalization
#' dd_loess = NormalizeBeta(dd, samples=c("dmso", "plx"), method="loess")
#' head(dd_loess)
#'
#'
#' @export

#===normalize function=====================================
NormalizeBeta <- function(beta, id = 1, method="cell_cycle",
                          posControl=NULL, samples=NULL){
  normalized = beta[, colnames(beta)[setdiff(1:ncol(beta), id)]]
  if(id==0) ids = rownames(beta) else ids = as.character(beta[,id])
  if(!is.null(samples)) normalized = normalized[, samples]
  normalized = as.matrix(normalized)
  if(method=="cell_cycle"){
    if(is.null(posControl)){
      Zuber_Essential = NULL
      data(Zuber_Essential)
      posControl=Zuber_Essential$GeneSymbol
    }
    idx = which(ids %in% toupper(posControl))
    if(length(idx)>0){
      mid = apply(normalized[idx,], 2, median, na.rm = TRUE)
      # mad = apply(normalized[idx,], 2, mad, na.rm = TRUE)
      mid = abs(mid - 0.1)
      normalized = t(t(normalized) / mid)
    }else{
      warning("No positive control genes are mapped !!!", call. = FALSE)
    }
  }
  if(method=="loess"){
    normalized = normalize.loess(normalized, log.it = FALSE, verbose=FALSE)
  }
  beta[, samples] = normalized
  return(beta)
}


#' normalize.loess
#'
#' Loess normalization method.
#'
#' @docType methods
#' @name normalize.loess
#' @rdname normalize.loess
#' @aliases loess.normalize
#'
#' @param mat A matrix with columns containing the values of the chips to normalize.
#' @param subset A subset of the data to fit a loess to.
#' @param epsilon A tolerance value (supposed to be a small value - used as a stopping criterion).
#' @param maxit Maximum number of iterations.
#' @param log.it Logical. If \code{TRUE} it takes the log2 of \code{mat}.
#' @param verbose Logical. If \code{TRUE} displays current pair of chip being worked on.
#' @param span Parameter to be passed the function \code{\link[stats]{loess}}
#' @param family.loess Parameter to be passed the function \code{\link[stats]{loess}}.
#' \code{"gaussian"} or \code{"symmetric"} are acceptable values for this parameter.
#' @param ... Any of the options of normalize.loess you would like to modify (described above).
#'
#' @return A matrix similar as \code{mat}.
#'
#' @author Wubing Zhang
#'
#' @seealso \code{\link{loess}}
#' @seealso \code{\link{NormalizeBeta}}
#'
#' @examples
#' beta = ReadBeta(mle.gene_summary)
#' beta_loess = normalize.loess(beta[,c("dmso", "plx")])
#'
#' @export
#'
normalize.loess <- function(mat, subset=sample(1:(dim(mat)[1]), min(c(5000, nrow(mat)))),
                            epsilon=10^-2, maxit=1, log.it=FALSE, verbose=TRUE, span=2/3,
                            family.loess="symmetric", ...){

  J <- dim(mat)[2]
  II <- dim(mat)[1]
  if(log.it){
    mat <- log2(mat)
  }

  change <- epsilon +1
  iter <- 0
  w <- c(0, rep(1,length(subset)), 0) ##this way we give 0 weight to the
  ##extremes added so that we can interpolate

  while(iter < maxit){
    iter <- iter + 1
    means <- matrix(0,II,J) ##contains temp of what we substract

    for (j in 1:(J-1)){
      for (k in (j+1):J){
        y <- mat[,j] - mat[,k]
        x <- (mat[,j] + mat[,k]) / 2
        index <- c(order(x)[1], subset, order(-x)[1])
        ##put endpoints in so we can interpolate
        xx <- x[index]
        yy <- y[index]
        aux <-loess(yy~xx, span=span, degree=1, weights=w, family=family.loess)
        aux <- predict(aux, data.frame(xx=x)) / J
        means[, j] <- means[, j] + aux
        means[, k] <- means[, k] - aux
        if (verbose)
          cat("Done with",j,"vs",k,"in iteration",iter,"\n")
      }
    }
    mat <- mat - means
    change <- max(colMeans((means[subset,])^2))

    if(verbose)
      cat(iter, change,"\n")

  }

  if ((change > epsilon) & (maxit > 1))
    warning(paste("No convergence after", maxit, "iterations.\n"))

  if(log.it) {
    return(2^mat)
  } else
    return(mat)
}
