#maplot

MA.plot <- function(beta, ctrlname="Control",treatname="Treatment", main="Negative control normalized", subset = sample(1:length(M), min(c(10000, length(M)))),
                    show.statistics = TRUE, span = 2/3, family.loess = "gaussian",
                    cex = 1, cex.lab=1.2, cex.axis=1, cex.main=1.2, plot.method = c("normal","smoothScatter","add"),
                    add.loess = TRUE, lwd = 1, lty = 1, loess.col = "red",filename=NULL){

  dd=beta
  loginfo(paste("MAplot for", main, "beta scores ..."))
  A = rowMeans(dd[,c(ctrlname, treatname)])
  M = rowMeans(dd[,treatname,drop=F])-rowMeans(dd[,ctrlname,drop=F])
  par(cex.axis=cex.axis, cex.lab=cex.lab,cex.main=cex.main)
  ma.plot(A, M, main=main, subset=subset, show.statistics=show.statistics, span=span,
          family.loess=family.loess, cex=cex, plot.method=plot.method, add.loess=add.loess,
          lwd=lwd, lty=lty, loess.col=loess.col)
  if(!is.null(filename)){
    png(filename,units = "in",width=500/100,height =450/100, res=300)
    par(cex.axis=cex.axis, cex.lab=cex.lab,cex.main=cex.main)
    ma.plot(A, M, main=main, subset=subset, show.statistics=show.statistics, span=span,
            family.loess=family.loess, cex=cex, plot.method=plot.method, add.loess=add.loess,
            lwd=lwd, lty=lty, loess.col=loess.col)
    dev.off()
  }
}
