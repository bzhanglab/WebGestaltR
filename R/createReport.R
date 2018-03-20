createReport <- function(hostName,outputDirectory,organism="hsapiens",timeStamp,enrichMethod="SEA",geneSet,geneSetDes,geneSetDAG,geneSetNet,interestingGeneMap,referenceGeneList,enrichedSig,enrichDatabase="geneontology_Biological_Process",enrichDatabaseFile=NULL,enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL,interestGene=NULL,interestGeneType=NULL,collapseMethod="mean",referenceGeneFile=NULL,referenceGene=NULL,referenceGeneType=NULL,referenceSet=NULL,minNum=10,maxNum=500,fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,topThr=10,dNum=20,perNum=1000,lNum=20,dagColor="binary"){

	outputHtmlFile <- file.path(outputDirectory,paste("Project_",timeStamp,sep=""),paste("Report_",timeStamp,".html",sep=""))

	# if hostname starts with "file://", it is used as WebGestaltReporter
	# all web assets are avaialble inside a parent directory called "assets"
	if(length(grep("file://", hostName, fixed=TRUE))==1){
		#file.symlink("../assets", file.path(outputDirectory, paste("Project_",timeStamp,sep=""),"assets"))
		#hostName <- "assets"
		hostName <- "https://s3-us-west-2.amazonaws.com/webgestalt/assets"
	}

	#htmlTitle(outputHtmlFile,hostName,organism,geneSetNet,geneSetDAG)
	if(organism!="others"){
		#createHTMLTab(outputHtmlFile,enrichedSig,interestingGeneMap)
		#####Summary Tab########
		tabsContent <- summaryTab_description(timeStamp,organism,interestGeneFile,interestGene,interestGeneType,enrichMethod,enrichDatabase,enrichDatabaseFile,enrichDatabaseType,enrichDatabaseDescriptionFile,interestingGeneMap,referenceGeneList,referenceGeneFile,referenceGene,referenceGeneType,referenceSet,minNum,maxNum,sigMethod,fdrThr,topThr,fdrMethod,enrichedSig,dNum,perNum,lNum,geneSet)

		###############Plot the mapping table######################
		tabsContent <- paste(tabsContent, mappingTableTab(interestingGeneMap), sep='\n')

		standardId <- interestingGeneMap$standardId
		###########GOSlim summary#########################
		if(standardId=="entrezgene"){
			tabsContent <- paste(tabsContent, goslimReportTab(timeStamp), sep='\n')
		}

		############Enrichment result##################
		if(!is.null(enrichedSig)){
			tabsContent <- paste(tabsContent, enrichResultTab_categoryViz(outputHtmlFile,organism,enrichMethod,fdrMethod,enrichedSig,dNum,geneSetDAG,geneSetDes,geneSetNet,outputDirectory,timeStamp,dagColor,hostName,interestingGeneMap,enrichDatabase), sep='\n')
		}

		template <- readLines(system.file("inst/templates/tab.mustache", package="WebGestaltR"))
		data <- list(hasEnrichedSig=!is.null(enrichedSig), idIsEntrezGene=interestingGeneMap$standardId=="entrezgene", tabsContent=tabsContent)
		bodyContent <- whisker.render(template, data)
	}else{
		###########Organism is others. No mapping information#############
		#############summary for the analysis###################
		bodyContent <- summaryTab_description(timeStamp,organism,interestGeneFile,interestGene,interestGeneType,enrichMethod,enrichDatabase,enrichDatabaseFile,enrichDatabaseType,enrichDatabaseDescriptionFile,interestingGeneMap,referenceGeneList,referenceGeneFile,referenceGene,referenceGeneType,referenceSet,minNum,maxNum,sigMethod,fdrThr,topThr,fdrMethod,enrichedSig,dNum,perNum,lNum,geneSet)

		##############Enrich Result################
		if(!is.null(enrichedSig)){
			bodyContent <- paste(bodyContent, enrichResult_Others(enrichMethod,enrichedSig,geneSetDes,fdrMethod,dNum), sep='\n')
		}
	}

	template <- readLines(system.file("inst/templates/template.mustache", package="WebGestaltR"))
	data <- list(hostName=hostName, organismIsOthers=organism=="others", geneSetNet=geneSetNet, geneSetDAG=geneSetDAG, bodyContent=bodyContent)
	cat(whisker.render(template, data), file=outputHtmlFile)
}
