#===normalize function=====================================
NormalizeBeta <- function(beta, method="cell_cycle", posControl=NULL, minus=0.6){
  loginfo("Normalize beta scores ...")
  if(method=="cell_cycle"){
    if(!is.null(posControl) && class(posControl)=="character" && file.exists(posControl)[1]){
      tmp = read.table(posControl, sep = "\t", header = FALSE)
      posControl = as.character(unlist(tmp))
    }else{
      posControl=essential_list
    }
    normalized = beta
    idx = which(normalized$Gene %in% posControl)
    mid = apply(normalized[idx, 2:(ncol(normalized)-1)],2,median)
    mid = abs(mid - minus)
    normalized[, 2:(ncol(normalized)-1)] = t(t(normalized[,2:(ncol(normalized)-1)]) / mid)
  }
  if(method=="loess"){
    normalized = beta
    normalized[,2:(ncol(normalized)-1)]=normalize.loess(normalized[,2:(ncol(normalized)-1)],log.it = FALSE)
  }
  return(normalized)
}
