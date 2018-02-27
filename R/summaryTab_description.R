summaryTab_description <- function(outputHtmlFile,timeStamp,organism,interestGeneFile,interestGene,interestGeneType,enrichMethod,enrichDatabase,enrichDatabaseFile,enrichDatabaseType,enrichDatabaseDescriptionFile,interestingGeneMap,referenceGene_List,referenceGeneFile,referenceGene,referenceGeneType,referenceSet,minNum,maxNum,sigMethod,fdrThr,topThr,fdrMethod,enrichedSig,dNum,perNum,lNum,geneSet){
	cat('<div id="summary">\n',file=outputHtmlFile,append=TRUE)
	#############summary for the analysis###################
	cat('<h3>Summary&nbsp&nbsp&nbsp(<a href="Project_',timeStamp,'.tar.gz" target="_blank" style="color:red"><u>Result Download</u></a>)</h3>\n',file=outputHtmlFile,append=TRUE,sep="")
	cat("<b>Enrich method:</b> ",enrichMethod,"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")

	cat("<b>Organism: </b>",organism,"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")

	###Introduce enrichdatabase
	if(organism!="others"){
		if(enrichDatabase=="others"){
			if(is.null(enrichDatabaseDescriptionFile)){
				cat("<b>Enrichment Categories: </b>",enrichDatabaseFile," <b>ID Type:</b>",enrichDatabaseType,"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
			}else{
				cat("<b>Enrichment Categories: </b>",enrichDatabaseFile," <b>ID Type:</b>",enrichDatabaseType," <b>Description File:</b>",enrichDatabaseDescriptionFile,"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
			}
		}else{
			cat("<b>Enrichment Categories: </b>",enrichDatabase,"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
		}
	}

	###Introduce interesting gene list
	if(organism!="others"){
		if(!is.null(interestGeneFile)){
			cat("<b>Interesting list: </b>",basename(interestGeneFile),". <b>ID type: </b>",interestGeneType,"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
		}else{
			cat("<b>Interesting list: </b> a R object. <b> ID type: </b>",interestGeneType,"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
		}

		standardId <- interestingGeneMap$standardId
		cat("The interesting list contains <b>",nrow(interestingGeneMap$mapped)+length(interestingGeneMap$unmapped),"</b> user IDs in which <b>",nrow(interestingGeneMap$mapped),"</b> user IDs are unambiguously mapped to <b>",length(unique(interestingGeneMap$mapped[,standardId])),"</b> unique ",standardId," IDs and <b>",length(interestingGeneMap$unmapped),"</b> user IDs can not be mapped to any ",standardId," ID.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")

		if(standardId=="entrezgene"){
			cat("The GO Slim summary are based upon the <b>",length(unique(interestingGeneMap$mapped[,standardId])),"</b> unique ",standardId," IDs.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
		}
	}else{
		if(!is.null(interestGeneFile)){
			cat("<b>Interesting list: </b>",basename(interestGeneFile),".<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
		}else{
			cat("<b>Interesting list: </b> a R object.<b><br/>\n",file=outputHtmlFile,append=TRUE,sep="")
		}
	}

	####Hard Code#####
	if(enrichMethod=="ORA"){
		sepcificParameterSummary_ORA(organism,outputHtmlFile,referenceGene_List,geneSet,referenceGeneFile,referenceGene,referenceGeneType,referenceSet,minNum,maxNum,sigMethod,fdrThr,topThr,fdrMethod,enrichedSig,dNum,interestingGeneMap)
	}

	if(enrichMethod=="GSEA"){
		sepcificParameterSummary_GSEA(organism,outputHtmlFile,interestingGeneMap,geneSet,minNum,maxNum,sigMethod,fdrThr,topThr,perNum,lNum,enrichedSig,dNum)
	}
	cat('</div>\n',file=outputHtmlFile,append=TRUE)  ##Summary END
}


sepcificParameterSummary_ORA <- function(organism,outputHtmlFile,referenceGene_List,geneSet,referenceGeneFile,referenceGene,referenceGeneType,referenceSet,minNum,maxNum,sigMethod,fdrThr,topThr,fdrMethod,enrichedSig,dNum,interestingGeneMap){
	###For others organisms, standardId is NULL
	if(organism!="others"){
		standardId <- interestingGeneMap$standardId
		interestGene_List <- unique(interestingGeneMap$mapped[,standardId])
		cat("Among <b>",length(interestGene_List),"</b> unique ",standardId," IDs, <b>", length(intersect(interestGene_List,intersect(referenceGene_List,geneSet[,3]))),"</b> IDs are annotated to the selected functional categories and also in the reference list, which are used for the enrichment analysis.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")

		##Introduce reference gene list
		if(!is.null(referenceGeneFile) || !is.null(referenceGene)){
			if(!is.null(referenceGeneFile)){
				cat("<b>Reference list: </b>",referenceGeneFile," <b>ID type: </b>",referenceGeneType,"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
			}else{
				cat("<b>Reference list: </b>a R object. <b>ID type: </b>",referenceGeneType,"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
			}
		}else{
			cat("<b>Reference list: </b> all mapped ",standardId," IDs from the selected platform ",referenceSet,"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
		}
		cat("The reference list can be mapped to <b>",length(referenceGene_List),"</b> ",standardId," IDs and <b> ",length(intersect(referenceGene_List,geneSet[,3])),"</b> IDs are annotated to the selected functional categories that are used as the reference for the enrichment analysis.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
	}else{  ###for others
		interestGene_List <- unique(interestingGeneMap[,1])
		cat("The file contains <b>",length(unique(interestGene_List)),"</b> user IDs (no ID mapping). All these IDs are used to perform the enrichment analysis.<br/>\n ",file=outputHtmlFile,append=TRUE,sep="")
		if(!is.null(referenceGeneFile)){
			cat("<b>Reference List: </b>",referenceGeneFile,"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
		}else{
			cat("<b>Reference List: </b>a R object<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
		}
		cat("<b>Total number of reference IDs: </b>",length(referenceGene_List),"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
	}

	cat("<b>Parameters for the enrichment analysis:</b>\n",file=outputHtmlFile,append=TRUE,sep="")
	cat("<ul><li><b>Minimum number of IDs in the category:</b>",minNum,"</li>\n",file=outputHtmlFile,append=TRUE,sep="")
	cat("<li><b>Maximum number of IDs in the category:</b>",maxNum,"</li>\n",file=outputHtmlFile,append=TRUE,sep="")
	cat("<li><b>FDR Method:</b>",fdrMethod,"</li>\n",file=outputHtmlFile,append=TRUE,sep="")
	if(sigMethod=="fdr"){
		cat("<li><b>Significance Level:</b> FDR<",fdrThr,"</li>\n",file=outputHtmlFile,append=TRUE,sep="")
	}

	if(sigMethod=="top"){
		cat("<li><b>Significance Level:</b> Top",topThr,"</li>\n",file=outputHtmlFile,append=TRUE,sep="")
	}
	cat("</ul>\n",file=outputHtmlFile,append=TRUE,sep="")
	if(!is.null(enrichedSig)){
		if(dNum>=nrow(enrichedSig)){
			cat("Based on the above parameters, <b>",nrow(enrichedSig),"</b> categories are identified as enriched categories and all are shown in this report.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
		}else{
			cat("Based on the above parameters, <b>",nrow(enrichedSig),"</b> categories are identified as enriched categories, in which <b>",dNum,"</b> most significant categories are shown in this report. All significant categories can be downloaded from the 'Result Download' link.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
		}
	}else{
		cat("Based on the above parameters, <b>No</b> category is identified as enriched category.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
	}
}


sepcificParameterSummary_GSEA <- function(organism,outputHtmlFile,interestingGeneMap,geneSet,minNum,maxNum,sigMethod,fdrThr,topThr,perNum,lNum,enrichedSig,dNum){
	if(organism!="others"){
		standardId <- interestingGeneMap$standardId
		interestGene_List <- unique(interestingGeneMap$mapped[,standardId])
		cat("Among the <b>",nrow(interestGene_List),"</b> unique ",standardId," IDs, <b>", length(intersect(interestGene_List,geneSet[,3])),"</b> IDs are annotated to the selected functional categories, which are used for the enrichment analysis.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
	}else{
		interestGene_List <- unique(interestingGeneMap[,1])
		cat("The file contains <b>",length(interestGene_List),"</b> user IDs (no ID mapping). All these IDs are used to perform the enrichment analysis.<br/>\n ",file=outputHtmlFile,append=TRUE,sep="")
	}

	cat("<b>Parameters for the enrichment analysis:</b>\n",file=outputHtmlFile,append=TRUE,sep="")
	cat("<ul><li><b>Minimum number of IDs in the category:</b>",minNum,"</li>\n",file=outputHtmlFile,append=TRUE,sep="")
	cat("<li><b>Maximum number of IDs in the category:</b>",maxNum,"</li>\n",file=outputHtmlFile,append=TRUE,sep="")
	if(sigMethod=="fdr"){
		cat("<li><b>Significance Level:</b> FDR<",fdrThr,"</li>\n",file=outputHtmlFile,append=TRUE,sep="")
	}
	if(sigMethod=="top"){
		cat("<li><b>Significance Level:</b> Top",topThr,"</li>\n",file=outputHtmlFile,append=TRUE,sep="")
	}
	cat("<li><b>Number of permutation:</b>",perNum,"</li>\n",file=outputHtmlFile,append=TRUE,sep="")
	cat("<li><b>Number of categories with the outputted leading edge IDs:</b>",lNum,"</li></ul>\n",file=outputHtmlFile,append=TRUE,sep="")

	if(!is.null(enrichedSig)){
		x <- enrichedSig[enrichedSig[,"NES"]>0,]
		y <- enrichedSig[enrichedSig[,"NES"]<0,]

		if(nrow(x)>0){
			if(dNum>=nrow(x)){
				cat("Based on the above parameters, <b>",nrow(x)," positive related </b>categories are identified as enriched categories and all are shown in this report.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
			}else{
				cat("Based on the above parameters, <b>",nrow(x)," positive related </b>categories are identified as enriched categories, in which <b>",dNum,"</b> most significant categories are shown in this report. All positive related significant categories can be downloaded from the 'Result Download' link.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
			}
		}else{
			cat("Based on the above parameters, <b>No positive related </b>category is identified as enriched category.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
		}

		if(nrow(y)>0){
			if(dNum>=nrow(y)){
				cat("Based on the above parameters, <b>",nrow(y)," negative related </b>categories are identified as enriched categories and all are shown in this report.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
			}else{
				cat("Based on the above parameters, <b>",nrow(y)," negative related </b>categories are identified as enriched categories, in which <b>",dNum,"</b> most significant categories are shown in this report. All positive related significant categories can be downloaded from the 'Result Download' link.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
			}
		}else{
			cat("Based on the above parameters, <b>No negative related </b>category is identified as enriched category.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
		}
	}else{
		cat("Based on the above parameters, <b>No</b> category is identified as enriched category.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
	}
}
