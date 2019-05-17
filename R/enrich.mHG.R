enrich.mHG = function(geneList, keytype = "Entrez",
                      type = "CORUM+GOBP+GOMF+GOCC+KEGG",
                      organism = 'hsa', pvalueCutoff = 0.05,
                      limit = c(3, 80),  gmtpath = NA){
  requireNamespace("clusterProfiler", quietly=TRUE) || stop("need clusterProfiler package")
  requireNamespace("data.table", quietly=TRUE) || stop("need data.table package")

  ## Prepare gene set annotation
  if(is.na(gmtpath)){
    msigdb = file.path(system.file("extdata", package = "MAGeCKFlute"),
                       paste0(organism, "_msig_entrez.gmt.gz"))
    gmtpath = gzfile(msigdb)
  }
  gene2path = ReadGMT(gmtpath, limit = limit)
  close(gmtpath)
  names(gene2path) = c("Gene","PathwayID", "PathwayName")
  gene2path$PathwayName = paste0(toupper(substr(gene2path$PathwayName, 0, 1)),
                                 tolower(substr(gene2path$PathwayName, 2,
                                                nchar(gene2path$PathwayName))))
  if(type != "All"){
    type = unlist(strsplit(type, "\\+"))
    idx = toupper(gsub("_.*", "", gene2path$PathwayID)) %in% toupper(type)
    gene2path = gene2path[idx, ]
  }
  gene2path = gene2path[!is.na(gene2path$Gene), ]
  idx = duplicated(gene2path$PathwayID)
  pathways = data.frame(PathwayID = gene2path$PathwayID[!idx],
                        PathwayName = gene2path$PathwayName[!idx],
                        stringsAsFactors = FALSE)

  ## Gene ID conversion
  if(keytype != "Entrez"){
    allsymbol = names(geneList)
    gene = TransGeneID(allsymbol, keytype, "Entrez", organism = organism)
    geneList = geneList[!duplicated(gene)]
    names(geneList) = gene[!duplicated(gene)]
  }
  gene = names(geneList)
  allsymbol = TransGeneID(gene, "Entrez", "Symbol", organism = organism)
  gene2path$Gene = as.character(gene2path$Gene)

  sMHG <- function(U,K,IDX){
    nseq <- which(IDX==1)
    k=sum(IDX)
    p_val <- sapply(nseq,function(x){
      subidx<-IDX[1:x]
      k_n <- sum(subidx)
      p<-1-phyper(q=k_n-1,k=x,m=K,n=U-K)
    })
    n_r <- nseq[which.min(p_val)]
    smHG <- min(p_val)
    return(smHG)
  }
  mHG <- function(W,K,U,IDX){
    N1=t(matrix(rep(seq(1,W),K),nrow = W))
    N2=matrix(rep(seq(1,K),W),nrow = K)
    N=U-(N1+N2)+1
    dk=(K-N2+1)/N
    dw=(W-N1+1)/N
    dk0=(K-seq(1,K)+1)/(U-seq(1,K)+1)
    dw0=(W-seq(1,W)+1)/(U-seq(1,W)+1)
    smhg =sMHG(U,K,IDX)
    mvalue <- function(Kx,Wy,smhg){
      if((1-phyper(q=Kx-1,m=K,k=Kx+Wy,n=W))<=smhg) return(0)
      if(Kx==0&&Wy==0) m=1
      if(Kx==0&&Wy>0) m=mvalue(Kx,Wy-1,smhg)*dw0[Wy]
      if(Wy==0&&Kx>0) m=mvalue(Kx-1,Wy,smhg)*dk0[Kx]
      if(Kx>0&&Wy>0) m=mvalue(Kx-1,Wy,smhg)*dk[Kx,Wy]+mvalue(Kx,Wy-1,smhg)*dw[Kx,Wy]
      return(m)
    }
    pval=1-mvalue(Kx=K,Wy=W,smhg)
    return(pval)
  }


  ## Start to do the hypergeometric test ##
  mHG <- function(pid){
    pGene = gene2path$Gene[gene2path$PathwayID==pid]
    IDX=as.numeric(pGene%in%gene)
    U=length(pGene)
    K=sum(IDX)
    W=U-K
    geneID = paste(pGene[which(idx==1)], collapse = "/")
    retr <- list(ID = NA, Description = NA, NES = NA,
                 pvalue = NA, GeneRatio = NA, BgRatio = NA,
                 geneID = NA, geneName = NA, Count = NA)
    if(k>=limit[1] & k<=limit[2] & q>2){
      pvalue = phyper(W,K,U,IDX)
      retr <- list(ID = pid, Description = pathways$PathwayName[pathways$PathwayID==pid],
                   NES = mean(geneList[pGene[which(idx==1)]]) * K^0.6,
                   pvalue = pvalue, GeneRatio = paste(q, sum(idx1), sep="/"),
                   BgRatio = paste(sum(idx1), length(pGene), sep="/"),
                   geneID = geneID, geneName = paste(allsymbol[pGene[which(idx==1)]],
                                                     collapse = "/"), Count = q)
    }
    return(retr)
  }

  ## Test using above function ##
  len = length(unique(intersect(gene, gene2path$Gene)))
  message("\t", len, " genes are mapped ...")
  res = sapply(pathways$PathwayID, mHG)
  res = as.data.frame(t(res), stringsAsFactors = FALSE)
  res = res[!is.na(res$ID), ]
  if(nrow(res)>0){
    res[, c(1:2, 5:8)] = matrix(unlist(res[, c(1:2, 5:8)]), ncol = 6)
    res[, c(3:4, 9)] = matrix(unlist(res[, c(3:4, 9)]), ncol = 3)
    res$p.adjust = p.adjust(res$pvalue, "BH")
    res$nLogpvalue = -log10(res$p.adjust)
    idx = which(res$p.adjust<=pvalueCutoff)
    if(length(idx)>0){
      res = res[idx, ]
      idx = c("ID", "Description", "NES", "pvalue", "p.adjust",
              "GeneRatio", "BgRatio", "geneID", "geneName", "Count")
      res = res[, idx]
    }else res=data.frame()
  }

  ## Create enrichResult object ##
  new("enrichResult",
      result         = res,
      pvalueCutoff   = pvalueCutoff,
      pAdjustMethod  = "BH",
      organism       = organism,
      ontology       = type,
      gene           = as.character(gene),
      keytype        = keytype)
}
