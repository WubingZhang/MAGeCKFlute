
BestRandom=function(d){
  den=density(d)
  denM=max(den$y)
  offset=seq(-1,1,0.01)
  M=10000000000
  for(i in offset){
    r <- rnorm(length(d), mean=0, sd=sd(d)+i)
    if(any(is.na(r))){next}
    denR=density(r)
    denRM=max(denR$y)
    if(abs(denM-denRM)<M){
      myR=r
      M=abs(denM-denRM)
    }
  }
  return(myR)
}


