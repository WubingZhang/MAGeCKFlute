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
#' @param beta Data frame, in which rows are EntrezID, columns are samples.
#' @param samples Character vector, specifying the samples in \code{beta} to be normalized.
#' If NULL (default), normalize beta score of all samples in \code{beta}.
#' @param method Character, one of 'cell_cycle'(default) and 'loess'.
#' @param posControl A file path or a character vector, specifying a list of gene entrezids as positive
#' controls used for cell cycle normalization
#' @param minus Numeric, scale for cell cycle normalization. Between 0 and 1.
#'
#' @return A data frame with same format as input data \code{beta}.
#'
#' @details In CRISPR screens, cells treated with different conditions (e.g., with or without
#' drug) may have different proliferation rates. So we defined a list of core essential genes,
#' which is equally negatively selected between samples with different proliferation rate.
#' Normalization of gene beta scores is performed using these essential genes. \code{cell_cycle}
#' in MAGeCKFlute normalizes the beta scores of all genes based on the median beta score of essential genes.
#' After normalization, the beta scores are comparable across samples. \code{loess} is another
#' optional normalization method, which is used to normalize array data before.
#'
#' @author Wubing Zhang
#'
#' @examples
#' data(MLE_Data)
#' # Read beta score from gene summary table in MAGeCK MLE results
#' dd = ReadBeta(MLE_Data, organism="hsa")
#' tmp = TransGeneID(rownames(dd), "Symbol", "Entrez")
#' dd = dd[!(duplicated(tmp)|is.na(tmp)), ]
#' rownames(dd) = tmp[!(duplicated(tmp)|is.na(tmp))]
#' samples=c("D7_R1", "D7_R2", "PLX7_R1", "PLX7_R2")
#' #Cell Cycle normalization
#' dd_essential = NormalizeBeta(dd, samples=samples, method="cell_cycle")
#' head(dd_essential)
#'
#' #Optional loess normalization
#' dd_loess = NormalizeBeta(dd, samples=samples, method="loess")
#' head(dd_loess)
#'
#'
#' @export

#===normalize function=====================================
NormalizeBeta <- function(beta, samples=NULL, method="cell_cycle", posControl=NULL, minus=0.2){
  message(Sys.time(), " # Normalize beta scores ...")
  if(is.null(samples)) samples = setdiff(colnames(beta))

  if(method=="cell_cycle"){
    if(!is.null(posControl) && class(posControl)=="character" && file.exists(posControl)[1]){
      tmp = read.table(posControl, sep = "\t", header = FALSE)
      posControl = as.character(unlist(tmp))
    }else{
      data(Zuber_Essential)
      posControl=Zuber_Essential
    }
    idx = which(rownames(beta) %in% posControl$EntrezID)
    normalized = as.matrix(beta[,samples])
    mid = apply(normalized[idx,], 2, median)
    mid = abs(mid - minus)
    normalized = t(t(normalized) / mid)
  }
  if(method=="loess"){
    normalized = as.matrix(beta[,samples])
    normalized = normalize.loess(normalized, log.it = FALSE, verbose=FALSE)
  }
  beta[,samples] = normalized

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
#' beta = ReadBeta(MLE_Data, organism="hsa")
#' beta_loess = normalize.loess(beta[,c("D7_R1", "D7_R2", "PLX7_R1", "PLX7_R2")])
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
