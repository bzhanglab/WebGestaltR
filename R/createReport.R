createReport <- function(hostName,outputDirectory,organism="hsapiens",timeStamp,enrichMethod="SEA",geneSet,geneSetDes,geneSetDAG,geneSetNet,interestingGeneMap,referenceGene_List,enrichedSig,enrichDatabase="geneontology_Biological_Process",enrichDatabaseFile=NULL,enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL,interestGene=NULL,interestGeneType=NULL,collapseMethod="mean",referenceGeneFile=NULL,referenceGene=NULL,referenceGeneType=NULL,referenceSet=NULL,minNum=10,maxNum=500,fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,topThr=10,dNum=20,perNum=1000,lNum=20,dagColor="binary"){
	 
	 outputHtmlFile <- file.path(outputDirectory,paste("Project_",timeStamp,sep=""),paste("Report_",timeStamp,".html",sep=""))
	 
   # if hostname starts with "file://", it is used as WebGestaltReporter
   # all web assets are avaialble inside a parent directory called "assets" 
   if(length(grep("file://", hostName, fixed=TRUE))==1){
     # create a symlink of assets in project directory
     # not sure if this will break on Windows
     file.symlink("../assets", file.path(outputDirectory, paste("Project_",timeStamp,sep=""),"assets"))
     hostName <- "assets"
   }

	 htmlTitle(outputHtmlFile,hostName,organism,geneSetNet,geneSetDAG)
	 if(organism!="others"){
			createHTMLTab(outputHtmlFile,enrichedSig,interestingGeneMap)
			#####Summary Tab########
			summaryTab_description(outputHtmlFile,timeStamp,organism,interestGeneFile,interestGene,interestGeneType,enrichMethod,enrichDatabase,enrichDatabaseFile,enrichDatabaseType,enrichDatabaseDescriptionFile,interestingGeneMap,referenceGene_List,referenceGeneFile,referenceGene,referenceGeneType,referenceSet,minNum,maxNum,sigMethod,fdrThr,topThr,fdrMethod,enrichedSig,dNum,perNum,lNum,geneSet)
		 
			###############Plot the mapping table######################
			mappingTableTab(outputHtmlFile,interestingGeneMap) 
			standardId <- interestingGeneMap$standardId
			###########GOSlim summary#########################
			if(standardId=="entrezgene"){
				goslimReportTab(outputHtmlFile,timeStamp)
			}
				
			############Enrichment result##################
			if(!is.null(enrichedSig)){
				cat('<div id="result">\n',file=outputHtmlFile,append=TRUE)
				dagInfo <- enrichResultTab_categoryViz(outputHtmlFile,organism,enrichMethod,fdrMethod,enrichedSig,dNum,geneSetDAG,geneSetDes,geneSetNet,outputDirectory,timeStamp,dagColor,hostName,interestingGeneMap,enrichDatabase)
				cat("</div>\n",file=outputHtmlFile,append=TRUE)
			}
			
			#######################################################
					
			cat("</div>\n",file=outputHtmlFile,append=TRUE)   ##This div is for the create tab function
			if(!is.null(enrichedSig) && !is.null(geneSetDAG)){	 		
				 jID <- paste("Project_",timeStamp,"_.json",sep="")
				 jF <- dagInfo$jF
				 style <- dagInfo$style
				 cat("<script>plotDAG('",jID,"','",jF,"','",style,"');</script>",file=outputHtmlFile,append=TRUE,sep="")
			}
			cat('<script>$( "#tabs" ).tabs();</script>',file=outputHtmlFile,append=TRUE)
				 		
		}else{
			###########Organism is others. No mapping information#############
			 #############summary for the analysis###################
	 		summaryTab_description(outputHtmlFile,timeStamp,organism,interestGeneFile,interestGene,interestGeneType,enrichMethod,enrichDatabase,enrichDatabaseFile,enrichDatabaseType,enrichDatabaseDescriptionFile,interestingGeneMap,referenceGene_List,referenceGeneFile,referenceGene,referenceGeneType,referenceSet,minNum,maxNum,sigMethod,fdrThr,topThr,fdrMethod,enrichedSig,dNum,perNum,lNum,geneSet)
			
			##############Enrich Result################
			if(!is.null(enrichedSig)){
				enrichResult_Others(outputHtmlFile,enrichMethod,enrichedSig,geneSetDes,fdrMethod,dNum)
			}
		}
	  cat('<iframe src="',file.path(hostName,"html","foot.html"),'" style="border:0px;height:130px;width:100%"></iframe>\n',file=outputHtmlFile,append=TRUE,sep="")
		cat("</body></html>",file=outputHtmlFile,append=TRUE)
}
