Cutoff_Calling=function(d,scale=F){
  param=1
  if(scale){param=round(length(d) / 16000,digits = 1)}
  Control_mean=0
  sorted_beta=sort(abs(d))
  temp=quantile(sorted_beta,0.68)
  temp_2=qnorm(0.84)
  cutoff=round(temp/temp_2,digits = 3)
  return(cutoff*param)
}
