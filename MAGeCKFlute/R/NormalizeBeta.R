#===normalize function=====================================
NormalizeBeta <- function(beta, method="cell_cycle", minus=0.6){
  data("Essential_list")
  loginfo("Normalize beta scores ...")
  if(method=="cell_cycle"){
    normalized = beta
    idx = which(normalized$Gene %in% essential_list)
    mid = apply(normalized[idx, 2:(ncol(normalized)-1)],2,median)
    mid = abs(mid - minus)
    normalized[, 2:(ncol(normalized)-1)] = t(t(normalized[,2:(ncol(normalized)-1)]) / mid)
  }
  if(method=="loess"){
    normalized = beta
    normalized[,2:(ncol(normalized)-1)]=normalize.loess(normalized[,2:(ncol(normalized)-1)],log.it = F)
  }
  return(normalized)
}
