#! /usr/bin/Rscript --vanilla
FluteMLE <- function(gene_summary, ctrlname="Control", treatname="Treatment",
                     prefix = "Test",scale_cutoff=1, adjust="none",
                     pvalueCutoff=0.05, enrich_kegg="ORT", organism="hsa",
                     gsea=F, loess=F, view_allpath=F, workspace="."){

	#=========Prepare the running environment=========
	{
	  loginfo("Create output dir and pdf file...")
	  out.dir=paste0(prefix, "_Flute_Results")
	  out.dir_sub=file.path(workspace,out.dir)
	  dir.create(file.path(out.dir_sub), showWarnings=FALSE)

	  output_pdf = paste0(prefix,"_Flute.mle_summary.pdf")
	  if(loess){
		pdf(file.path(out.dir_sub, output_pdf),width=15,height = 7)
	  }else{
		pdf(file.path(out.dir_sub, output_pdf),width=9,height = 6)
	  }
	}

	#=========Beta Score Preparation=========================
	{
	  dd = ReadBeta(gene_summary, ctrlName = ctrlname,
	                treatName = treatname, organism=organism)

	  if(loess){ dd_loess = NormalizeBeta(dd, method="loess") }
	  dd_essential = NormalizeBeta(dd, method="cell_cycle",minus=0.6)
	}


	# ========distribution of all genes================================
	{
	  outputDir1 <- file.path(out.dir_sub,"Distribution_of_BetaScores")
	  dir.create(outputDir1, showWarnings=FALSE)
	  #Negative normalization
	  P1=Violin.plot(dd, main="Negative control normalized",
	                 filename=file.path(outputDir1,
	                                    "violin_all_negative_normalized.png"))
	  P2=Density.plot(dd, main="Negative control normalized",
	                  filename=file.path(outputDir1,
	                                     "density_all_negative_normalized.png"))
	  P3=Density.diff(dd, main="Negative control normalized",
	                  filename=file.path(
	                    outputDir1, "density_all_treat-ctrl_negative_normalized.png"))
	  #Essential gene normalization
	  P4=Violin.plot(dd_essential, main="Cell cycle normalized",
	                 filename=file.path(outputDir1,
	                                    "violin_all_essential_normalized.png"))
	  P5=Density.plot(dd_essential, main="Cell cycle normalized",
	                  filename=file.path(
	                    outputDir1,"density_all_essential_normalized.png"))
	  P6=Density.diff(dd_essential, main="Cell cycle normalized",
	                  filename=file.path(
	                    outputDir1, "density_all_treat-ctrl_essential_normalized.png"))
	  #Loess normalization
	  if(loess){
		P7=Violin.plot(dd_loess,main="Loess normalized",
		               filename=file.path(
		                 outputDir1,"violin_all_loess_normalized.png"))
		P8=Density.plot(dd_loess,main="Loess normalized",
		                filename=file.path(
		                  outputDir1,"density_all_loess_normalized.png"))
		P9=Density.diff(dd_loess,main="Loess normalized",
		                filename=file.path(
		                  outputDir1,"density_all_treat-ctrl_loess_normalized.png"))
		grid.arrange(P1,P4,P7,P2,P5,P8, ncol = 3)
	  }else{grid.arrange(P1,P4,P2,P5, ncol = 2)}
	}

	#========MAplot of treatment and control beta scores==================
	{
	  outputDir2 <- file.path(out.dir_sub,"MAplot")
	  dir.create(outputDir2, showWarnings=FALSE)

	  if(loess){
		par(mfrow=c(2,3))
		plot.new()
		plot.new()
		vp.Top <- viewport(height=unit(0.5, "npc"), width=unit(1, "npc"),
		                   just=c("bottom"), y=0.5, x=0.5)
		plot.new()
		grid.arrange(P3,P6,P9,ncol = 3, vp=vp.Top,newpage = F)
	  }else{
		par(mfrow=c(2,2))
		plot.new()
		plot.new()
		vp.Top <- viewport(height=unit(0.5, "npc"), width=unit(1, "npc"),
		                   just=c("bottom"), y=0.5, x=0.5)
		grid.arrange(P3,P6,ncol = 2, vp=vp.Top,newpage = F)
	  }

	  MA.plot(dd, cex=1, main="Negative control normalized",
	          filename=file.path(outputDir2,"maplot_negative_normalized.png"))
	  MA.plot(dd_essential, cex=1, main="Cell cycle normalized",
	          filename=file.path(outputDir2,"maplot_essential_normalized.png"))
	  if(loess){MA.plot(dd_loess, cex=1, main="Loess normalized",
	                    filename=file.path(outputDir2,"maplot_loess_normalized.png"))}
	}


	#=============distribution of essential genes====================
	{
	  outputDir3=file.path(out.dir_sub,"Linear_Fitting_of_BetaScores")
	  dir.create(outputDir3, showWarnings=FALSE)

	  idx=which(dd$Gene %in% essential_list)
	  #Negative control normalized
	  P1=Violin.plot(dd[idx,],ylab="Essential.B.S.",
	                 main="Negative control normalized",
	                 filename=file.path(
	                   outputDir1,"violin_ess_negative_normalized.png"))
	  P2=Density.plot(dd[idx,],xlab="Essential.B.S.",
	                  main="Negative control normalized",
	                  filename=file.path(
	                    outputDir1,"density_ess_negative_normalized.png"))
	  P3=CellCycleFit(dd,ylab="Beta Score",
	                  main="Negative control normalized",
	                  filename=file.path(
	                    outputDir3,"Linear_all_negative_normalized.png"))
	  P4=CellCycleFit(dd[idx,],ylab="Essential.B.S.",
	                  main="Negative control normalized",
	                  filename=file.path(
	                    outputDir3,"Linear_ess_negative_normalized.png"))
	  #Essential normalized
	  P5=Violin.plot(dd_essential[idx,],ylab="Essential.B.S.",
	                 main="Cell cycle  normalized",
	                 filename=file.path(
	                   outputDir1,"violin_ess_essential_normalized.png"))
	  P6=Density.plot(dd_essential[idx,],xlab="Essential.B.S.",
	                  main="Cell cycle  normalized",
	                  filename=file.path(
	                    outputDir1,"density_ess_essential_normalized.png"))
	  P7=CellCycleFit(dd_essential,ylab="Beta Score",
	                  main="Cell cycle  normalized",
	                  filename=file.path(
	                    outputDir3,"Linear_all_essential_normalized.png"))
	  P8=CellCycleFit(dd_essential[idx,],ylab="Essential.B.S.",
	                  main="Cell cycle  normalized",
	                  filename=file.path(
	                    outputDir3,"Linear_ess_essential_normalized.png"))

	  #loess normalized
	  if(loess){
		P9=Violin.plot(dd_loess[idx,],ylab="Essential.B.S.",
		               main="Loess  normalized",
		               filename=file.path(
		                 outputDir1,"violin_ess_loess_normalized.png"))
		P10=Density.plot(dd_loess[idx,],xlab="Essential.B.S.",
		                 main="Loess  normalized",
		                 filename=file.path(
		                   outputDir1,"density_ess_loess_normalized.png"))
		P11=CellCycleFit(dd_loess,ylab="Beta Score",
		                 main="Loess  normalized",
		                 filename=file.path(
		                   outputDir3,"Linear_all_loess_normalized.png"))
		P12=CellCycleFit(dd_loess[idx,],ylab="Essential.B.S.",
		                 main="Loess  normalized",
		                 filename=file.path(
		                   outputDir3,"Linear_ess_loess_normalized.png"))

		grid.arrange(P1,P5,P9,P2,P6,P10, ncol=3)
		grid.arrange(P3,P7,P11,P4,P8,P12, ncol=3)
	  }else{
		grid.arrange(P1,P5,P2,P6, ncol=2)
		grid.arrange(P3,P7,P4,P8, ncol=2)
	  }
	}


	# =========drug-targeted genes=================================
	{
	  outputDir4 = file.path(out.dir_sub,"Scatter_Treat_Ctrl")
	  dir.create(outputDir4, showWarnings=FALSE)

	  #Negative normalized
	  ddAB=GroupAB(dd, scale_cutoff=scale_cutoff,
	               filename=file.path(
	                 outputDir4,"GroupAB_detail_negative_normalized.txt"))
	  P1=ScatterAB(ddAB,main="Negative control normalized",
	               scale_cutoff=scale_cutoff,
	               filename=file.path(
	                 outputDir4,"Scatter_Treat-Ctrl_negative_normalized.png"))
	  P2=RankAB(ddAB,genelist=c(),main="Negative control normalized",
	            scale_cutoff=scale_cutoff,
	            filename=file.path(
	              outputDir4,"Rank_Treat-Ctrl_negative_normalized.png"))
	  #Essential normalized
	  ddAB_essential=GroupAB(dd_essential, scale_cutoff=scale_cutoff,
	                         filename=file.path(
	                           outputDir4,"GroupAB_detail_essential_normalized.txt"))
	  P3=ScatterAB(ddAB_essential,main="Cell cycle normalized",
	               scale_cutoff=scale_cutoff,
	               filename=file.path(
	                 outputDir4,"Scatter_Treat-Ctrl_essential_normalized.png"))
	  P4=RankAB(ddAB_essential,genelist=c(),main="Cell cycle  normalized",
	            scale_cutoff=scale_cutoff,
	            filename=file.path(
	              outputDir4,"Rank_Treat-Ctrl_essential_normalized.png"))

	  #Loess normalized
	  if(loess){
		ddAB_loess=GroupAB(dd_loess, scale_cutoff=scale_cutoff,
		                   filename=file.path(
		                     outputDir4,"GroupAB_detail_loess_normalized.txt"))
		P5=ScatterAB(ddAB_loess,main="Loess normalized",
		             scale_cutoff=scale_cutoff,
		             filename=file.path(
		               outputDir4,"Scatter_Treat-Ctrl_loess_normalized.png"))
		P6=RankAB(ddAB_loess,genelist=c(),main="Loess  normalized",
		          scale_cutoff=scale_cutoff,
		          filename=file.path(
		            outputDir4,"Rank_Treat-Ctrl_loess_normalized.png"))

		grid.arrange(P1,P3,P5,P2,P4,P6, ncol=3)
	  }else{grid.arrange(P1,P3,P2,P4, ncol=2)}
	}

	#==========enrichment AB group genes=========================
	{
	  outputDir5 = file.path(out.dir_sub,"Enrichment_Treat-Ctrl/")
	  dir.create(outputDir5, showWarnings=FALSE)

	  E1 = EnrichAB(ddAB,pvalue=pvalueCutoff,enrich_method = enrich_kegg,
	                organism=organism,
	                adjust=adjust, gsea=gsea,
	                filename="Negative_ctrl_normalized", out.dir=outputDir5)
	  E2 = EnrichAB(ddAB_essential,pvalue=pvalueCutoff,
	                enrich_method = enrich_kegg, organism=organism,
	                adjust=adjust, gsea=gsea,
	                filename="Essential_normalized", out.dir=outputDir5)

	  if(loess){
  		E3 = EnrichAB(ddAB_loess,pvalue=pvalueCutoff,
  		              enrich_method = enrich_kegg, organism=organism,
  		              adjust=adjust, gsea=gsea,filename="Loess_normalized",
  		              out.dir=outputDir5)

  		grid.arrange(E1$keggA$gridPlot, E2$keggA$gridPlot, E3$keggA$gridPlot,
  		             E1$bpA$gridPlot, E2$bpA$gridPlot, E3$bpA$gridPlot,
  		             ncol = 3)
  		if(gsea) grid.arrange(E1$gseA$gridPlot, E2$gseA$gridPlot,
  		                      E3$gseA$gridPlot, E1$gseA$gseaplot,
  		                      E2$gseA$gseaplot, E3$gseA$gseaplot, ncol = 3)
  		grid.arrange(E1$keggB$gridPlot, E2$keggB$gridPlot, E3$keggB$gridPlot,
  		             E1$bpB$gridPlot, E2$bpB$gridPlot,
  		             E3$bpB$gridPlot, ncol = 3)
  		if(gsea)grid.arrange(E1$gseB$gridPlot, E2$gseB$gridPlot,
  		                     E3$gseB$gridPlot, E1$gseB$gseaplot,
  		                     E2$gseB$gseaplot, E3$gseB$gseaplot, ncol = 3)
	  }else{
      grid.arrange(E1$keggA$gridPlot, E2$keggA$gridPlot, E1$bpA$gridPlot,
                   E2$bpA$gridPlot, ncol = 2)
  		if(gsea) grid.arrange(E1$gseA$gridPlot, E2$gseA$gridPlot,
  		                      E1$gseA$gseaplot, E2$gseA$gseaplot, ncol = 2)
  		grid.arrange(E1$keggB$gridPlot, E2$keggB$gridPlot, E1$bpB$gridPlot,
  		             E2$bpB$gridPlot, ncol = 2)
  		if(gsea) grid.arrange(E1$gseB$gridPlot, E2$gseB$gridPlot,
  		                      E1$gseB$gseaplot, E2$gseB$gseaplot, ncol = 2)
	  }
	}

	#======Pathview for GroupA and GroupB significant enriched pathways
	{
	  dir.create(file.path(out.dir_sub,"Pathview_Treat_Ctrl"),
	             showWarnings=FALSE)
	  #Pathway view for top 4 pathway
	  pathPng2Pdf(dd, E1$keggA$enrichRes@result$ID, "Group A",
	              "Negative control normalized",organism=organism,
	              view_allpath=view_allpath, output=file.path(
	                out.dir_sub,"Pathview_Treat_Ctrl"))
	  pathPng2Pdf(dd, E1$keggB$enrichRes@result$ID, "Group B",
	              "Negative control normalized",organism=organism,
	              view_allpath=view_allpath, output=file.path(
	                out.dir_sub,"Pathview_Treat_Ctrl"))
	  pathPng2Pdf(dd_essential, E2$keggA$enrichRes@result$ID, "Group A",
	              "Cell cycle normalized",organism=organism,
	              view_allpath=view_allpath, output=file.path(
	                out.dir_sub,"Pathview_Treat_Ctrl"))
	  pathPng2Pdf(dd_essential, E2$keggB$enrichRes@result$ID, "Group B",
	              "Cell cycle normalized",organism=organism,
	              view_allpath=view_allpath, output=file.path(
	                out.dir_sub,"Pathview_Treat_Ctrl"))
	  if(loess){
		pathPng2Pdf(dd_loess, E3$keggA$enrichRes@result$ID, "Group A",
		            "Loess normalized",organism=organism,
		            view_allpath=view_allpath, output=file.path(
		              out.dir_sub,"Pathview_Treat_Ctrl"))
		pathPng2Pdf(dd_loess, E3$keggB$enrichRes@result$ID, "Group B",
		            "Loess normalized",organism=organism,
		            view_allpath=view_allpath, output=file.path(
		              out.dir_sub,"Pathview_Treat_Ctrl"))
	  }
	}

	# ===============9 squares=====================================
	{
	  dir.create(file.path(out.dir_sub,"Scatter_9Square"),
	             showWarnings=FALSE)
	  dd1=dd[,c("Gene","Treatment","Control","ENTREZID")]
	  P1=Square.plot(dd1,main="Negative control normalized",
	                 scale_cutoff=scale_cutoff, filename="negative_normalized",
	                 out.dir=file.path(out.dir_sub,"Scatter_9Square/"))

	  dd2=dd_essential[,c("Gene","Treatment","Control","ENTREZID")]
	  P2=Square.plot(dd2,main="Cell cycle normalized", scale_cutoff=scale_cutoff,
	                 filename="essential_normalized",
	                 out.dir=file.path(out.dir_sub,"Scatter_9Square/"))

	  if(loess){
  		dd3=dd_loess[,c("Gene","Treatment","Control","ENTREZID")]
  		P3=Square.plot(dd3,main="Loess normalized", scale_cutoff=scale_cutoff,
  		               filename="loess_normalized",
  		               out.dir=file.path(out.dir_sub,"Scatter_9Square/"))
  		grid.arrange(P1$p,P2$p,P3$p, ncol = 3)
	  }else{grid.arrange(P1$p,P2$p, ncol = 2, bottom="", left="",
	                     right="", top="")}
	}

	#==============9 Square grouped gene enrichment======================
	{
	  dir.create(file.path(out.dir_sub,"Enrichment_9Square"),
	             showWarnings=FALSE)
	  E1 = EnrichSquare(P1$dd1,pvalue=pvalueCutoff,adjust=adjust,
	                    enrich_method=enrich_kegg,organism=organism,
						 filename="negative_normalized", out.dir=
						   file.path(out.dir_sub,"Enrichment_9Square"))
	  E2 = EnrichSquare(P2$dd1,pvalue=pvalueCutoff,adjust=adjust,
	                    enrich_method=enrich_kegg,organism=organism,
						 filename="essential_normalized", out.dir=
						   file.path(out.dir_sub,"Enrichment_9Square"))

	  if(loess){
		E3 = EnrichSquare(P3$dd1,pvalue=pvalueCutoff,adjust=adjust,
		                  enrich_method=enrich_kegg, organism=organism,
						   filename="loess_normalized", out.dir=
						     file.path(out.dir_sub,"Enrichment_9Square"))

		grid.arrange(E1$kegg1$gridPlot, E2$kegg1$gridPlot, E3$kegg1$gridPlot,
		             E1$bp1$gridPlot, E2$bp1$gridPlot, E3$bp1$gridPlot,
		             ncol = 3)
		grid.arrange(E1$kegg2$gridPlot, E2$kegg2$gridPlot, E3$kegg2$gridPlot,
		             E1$bp2$gridPlot, E2$bp2$gridPlot, E3$bp2$gridPlot, ncol = 3)
		grid.arrange(E1$kegg3$gridPlot, E2$kegg3$gridPlot, E3$kegg3$gridPlot,
		             E1$bp3$gridPlot, E2$bp3$gridPlot, E3$bp3$gridPlot, ncol = 3)
		grid.arrange(E1$kegg4$gridPlot, E2$kegg4$gridPlot,
		             E3$kegg4$gridPlot, E1$bp4$gridPlot,
		             E2$bp4$gridPlot, E3$bp4$gridPlot, ncol = 3)
		grid.arrange(E1$kegg13$gridPlot, E2$kegg13$gridPlot,
		             E3$kegg13$gridPlot, E1$bp13$gridPlot,
		             E2$bp13$gridPlot, E3$bp13$gridPlot, ncol = 3)
		grid.arrange(E1$kegg14$gridPlot, E2$kegg14$gridPlot,
		             E3$kegg14$gridPlot, E1$bp14$gridPlot,
		             E2$bp14$gridPlot, E3$bp14$gridPlot, ncol = 3)
		grid.arrange(E1$kegg23$gridPlot, E2$kegg23$gridPlot,
		             E3$kegg23$gridPlot, E1$bp23$gridPlot,
		             E2$bp23$gridPlot, E3$bp23$gridPlot, ncol = 3)
		grid.arrange(E1$kegg24$gridPlot, E2$kegg24$gridPlot,
		             E3$kegg24$gridPlot, E1$bp24$gridPlot,
		             E2$bp24$gridPlot, E3$bp24$gridPlot, ncol = 3)
		grid.arrange(E1$kegg1234$gridPlot, E2$kegg1234$gridPlot,
		             E3$kegg1234$gridPlot, E1$bp1234$gridPlot,
		             E2$bp1234$gridPlot, E3$bp1234$gridPlot, ncol = 3)
	  }else{
		grid.arrange(E1$kegg1$gridPlot, E2$kegg1$gridPlot, E1$bp1$gridPlot,
		             E2$bp1$gridPlot, ncol = 2)
		grid.arrange(E1$kegg2$gridPlot, E2$kegg2$gridPlot, E1$bp2$gridPlot,
		             E2$bp2$gridPlot, ncol = 2)
		grid.arrange(E1$kegg3$gridPlot, E2$kegg3$gridPlot, E1$bp3$gridPlot,
		             E2$bp3$gridPlot, ncol = 2)
		grid.arrange(E1$kegg4$gridPlot, E2$kegg4$gridPlot, E1$bp4$gridPlot,
		             E2$bp4$gridPlot, ncol = 2)
		grid.arrange(E1$kegg13$gridPlot, E2$kegg13$gridPlot, E1$bp13$gridPlot,
		             E2$bp13$gridPlot, ncol = 2)
		grid.arrange(E1$kegg14$gridPlot, E2$kegg14$gridPlot, E1$bp14$gridPlot,
		             E2$bp14$gridPlot, ncol = 2)
		grid.arrange(E1$kegg23$gridPlot, E2$kegg23$gridPlot, E1$bp23$gridPlot,
		             E2$bp23$gridPlot, ncol = 2)
		grid.arrange(E1$kegg24$gridPlot, E2$kegg24$gridPlot, E1$bp24$gridPlot,
		             E2$bp24$gridPlot, ncol = 2)
		grid.arrange(E1$kegg1234$gridPlot, E2$kegg1234$gridPlot,
		             E1$bp1234$gridPlot, E2$bp1234$gridPlot, ncol = 2)
	  }
	}

	#========Pathway view for 9 square genesets================
	{
	  dir.create(file.path(out.dir_sub,"Pathview_9Square"),
	             showWarnings=FALSE)

	  pathPng2Pdf(dd, E1$kegg1$enrichRes@result$ID,"Group 1",
	              "Negative control normalized",organism=organism,
	              view_allpath=view_allpath, output=
	                file.path(out.dir_sub,"Pathview_9Square"))
	  pathPng2Pdf(dd_essential, E2$kegg1$enrichRes@result$ID,"Group 1",
	              "Cell cycle normalized",organism=organism,
	              view_allpath=view_allpath, output=
	                file.path(out.dir_sub,"Pathview_9Square"))

	  pathPng2Pdf(dd, E1$kegg2$enrichRes@result$ID,"Group 2",
	              "Negative control normalized",organism=organism,
	              view_allpath=view_allpath, output=
	                file.path(out.dir_sub,"Pathview_9Square"))
	  pathPng2Pdf(dd_essential, E2$kegg2$enrichRes@result$ID,"Group 2",
	              "Cell cycle normalized",organism=organism,
	              view_allpath=view_allpath, output=
	                file.path(out.dir_sub,"Pathview_9Square"))

	  pathPng2Pdf(dd, E1$kegg3$enrichRes@result$ID,"Group 3",
	              "Negative control normalized",organism=organism,
	              view_allpath=view_allpath, output=
	                file.path(out.dir_sub,"Pathview_9Square"))
	  pathPng2Pdf(dd_essential, E2$kegg3$enrichRes@result$ID,"Group 3",
	              "Cell cycle normalized",organism=organism,
	              view_allpath=view_allpath, output=
	                file.path(out.dir_sub,"Pathview_9Square"))

	  pathPng2Pdf(dd, E1$kegg4$enrichRes@result$ID,"Group 4",
	              "Negative control normalized",organism=organism,
	              view_allpath=view_allpath, output=
	                file.path(out.dir_sub,"Pathview_9Square"))
	  pathPng2Pdf(dd_essential, E2$kegg4$enrichRes@result$ID,"Group 4",
	              "Cell cycle normalized",organism=organism,
	              view_allpath=view_allpath, output=
	                file.path(out.dir_sub,"Pathview_9Square"))

	  pathPng2Pdf(dd, E1$kegg13$enrichRes@result$ID,"Group 1 & Group 3",
	              "Negative control normalized",organism=organism,
	              view_allpath=view_allpath, output=
	                file.path(out.dir_sub,"Pathview_9Square"))
	  pathPng2Pdf(dd_essential, E2$kegg13$enrichRes@result$ID,
	              "Group 1 & Group 3","Cell cycle normalized",
	              organism=organism, view_allpath=view_allpath,
	              output=file.path(out.dir_sub,"Pathview_9Square"))

	  pathPng2Pdf(dd, E1$kegg14$enrichRes@result$ID,"Group 1 & Group 4",
	              "Negative control normalized",organism=organism,
	              view_allpath=view_allpath, output=
	                file.path(out.dir_sub,"Pathview_9Square"))
	  pathPng2Pdf(dd_essential, E2$kegg14$enrichRes@result$ID,
	              "Group 1 & Group 4","Cell cycle normalized",
	              organism=organism, view_allpath=view_allpath,
	              output=file.path(out.dir_sub,"Pathview_9Square"))

	  pathPng2Pdf(dd, E1$kegg23$enrichRes@result$ID,"Group 2 & Group 3",
	              "Negative control normalized",organism=organism,
	              view_allpath=view_allpath, output=
	                file.path(out.dir_sub,"Pathview_9Square"))
	  pathPng2Pdf(dd_essential, E2$kegg23$enrichRes@result$ID,
	              "Group 2 & Group 3","Cell cycle normalized",
	              organism=organism, view_allpath=view_allpath,
	              output=file.path(out.dir_sub,"Pathview_9Square"))

	  pathPng2Pdf(dd, E1$kegg24$enrichRes@result$ID,
	              "Group 2 & Group 4","Negative control normalized",
	              organism=organism, view_allpath=view_allpath,
	              output=file.path(out.dir_sub,"Pathview_9Square"))
	  pathPng2Pdf(dd_essential, E2$kegg24$enrichRes@result$ID,
	              "Group 2 & Group 4","Cell cycle normalized",
	              organism=organism, view_allpath=view_allpath,
	              output=file.path(out.dir_sub,"Pathview_9Square"))

	  pathPng2Pdf(dd, E1$kegg1234$enrichRes@result$ID,
	              "Group 1 & Group 2 & Group 3 & Group 4",
	              "Negative control normalized",organism=organism,
	              view_allpath=view_allpath, output=
	                file.path(out.dir_sub,"Pathview_9Square"))
	  pathPng2Pdf(dd_essential, E2$kegg1234$enrichRes@result$ID,
	              "Group 1 & Group 2 & Group 3 & Group 4",
	              "Cell cycle normalized",organism=organism,
	              view_allpath=view_allpath, output=
	                file.path(out.dir_sub,"Pathview_9Square"))

	  if(loess){
		pathPng2Pdf(dd_loess, E3$kegg1$enrichRes@result$ID,"Group 1",
		            "Loess normalized",organism=organism,
		            view_allpath=view_allpath, output=
		              file.path(out.dir_sub,"Pathview_9Square"))
		pathPng2Pdf(dd_loess, E3$kegg2$enrichRes@result$ID,"Group 2",
		            "Loess normalized",organism=organism,
		            view_allpath=view_allpath, output=
		              file.path(out.dir_sub,"Pathview_9Square"))
		pathPng2Pdf(dd_loess, E3$kegg3$enrichRes@result$ID,"Group 3",
		            "Loess normalized",organism=organism,
		            view_allpath=view_allpath, output=
		              file.path(out.dir_sub,"Pathview_9Square"))
		pathPng2Pdf(dd_loess, E3$kegg4$enrichRes@result$ID,"Group 4",
		            "Loess normalized",organism=organism,
		            view_allpath=view_allpath, output=
		              file.path(out.dir_sub,"Pathview_9Square"))
		pathPng2Pdf(dd_loess, E3$kegg13$enrichRes@result$ID,
		            "Group 1 & Group 3","Loess normalized",organism=organism,
		            view_allpath=view_allpath, output=
		              file.path(out.dir_sub,"Pathview_9Square"))
		pathPng2Pdf(dd_loess, E3$kegg14$enrichRes@result$ID,
		            "Group 1 & Group 4","Loess normalized",organism=organism,
		            view_allpath=view_allpath, output=
		              file.path(out.dir_sub,"Pathview_9Square"))
		pathPng2Pdf(dd_loess, E3$kegg23$enrichRes@result$ID,
		            "Group 2 & Group 3","Loess normalized",organism=organism,
		            view_allpath=view_allpath, output=
		              file.path(out.dir_sub,"Pathview_9Square"))
		pathPng2Pdf(dd_loess, E3$kegg24$enrichRes@result$ID,
		            "Group 2 & Group 4","Loess normalized",organism=organism,
		            view_allpath=view_allpath, output=
		              file.path(out.dir_sub,"Pathview_9Square"))
		pathPng2Pdf(dd_loess, E3$kegg1234$enrichRes@result$ID,
		            "Group 1 & Group 2 & Group 3 & Group 4","Loess normalized",
		            organism=organism,view_allpath=view_allpath,
		            output=file.path(out.dir_sub,"Pathview_9Square"))
	  }
	}
	dev.off()
	#===============================================
}
