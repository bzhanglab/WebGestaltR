createReport <- function(hostName, outputDirectory, organism="hsapiens", timeStamp, enrichMethod, geneSet, geneSetDes, geneSetDag, geneSetNet, interestingGeneMap, referenceGeneList, enrichedSig, geneTables, background, enrichDatabase="geneontology_Biological_Process", enrichDatabaseFile=NULL, enrichDatabaseType=NULL, enrichDatabaseDescriptionFile=NULL, interestGeneFile=NULL, interestGene=NULL, interestGeneType=NULL, collapseMethod="mean", referenceGeneFile=NULL, referenceGene=NULL, referenceGeneType=NULL, referenceSet=NULL, minNum=10, maxNum=500, fdrMethod="BH", sigMethod="fdr", fdrThr=0.05, topThr=10, dNum=20, perNum=1000, lNum=20, dagColor="binary"){

	outputHtmlFile <- file.path(outputDirectory,paste("Project_",timeStamp,sep=""),paste("Report_",timeStamp,".html",sep=""))

	# if hostname starts with "file://", it is used as WebGestaltReporter
	# all web assets are avaialble inside a parent directory called "assets"
	## TODO: FIXME
	if(length(grep("file://", hostName, fixed=TRUE))==1){
		#file.symlink("../assets", file.path(outputDirectory, paste("Project_",timeStamp,sep=""),"assets"))
		#hostName <- "assets"
		hostName <- "https://s3-us-west-2.amazonaws.com/webgestalt/assets"
	}

	numAnnoRefUserId <- NULL
	dagJson <- NULL
	if(organism!="others"){
		#####Summary Tab########
		tabsContent <- summaryDescription(timeStamp,organism,interestGeneFile,interestGene,interestGeneType,enrichMethod,enrichDatabase,enrichDatabaseFile,enrichDatabaseType,enrichDatabaseDescriptionFile,interestingGeneMap,referenceGeneList,referenceGeneFile,referenceGene,referenceGeneType,referenceSet,minNum,maxNum,sigMethod,fdrThr,topThr,fdrMethod,enrichedSig,dNum,perNum,lNum,geneSet)

		###############Plot the mapping table######################
		# tabsContent <- paste(tabsContent, mappingTableTab(interestingGeneMap), sep='\n')

		standardId <- interestingGeneMap$standardId
		if (enrichMethod == 'ORA') {
			interestGeneList <- unique(interestingGeneMap$mapped[,standardId])
			numAnnoRefUserId <- length(intersect(interestGeneList,intersect(referenceGeneList,geneSet[,3])))
		}
		###########GOSlim summary#########################
		if(standardId=="entrezgene"){
			tabsContent <- paste(tabsContent, goSlimReport(timeStamp), sep='\n')
		}

		############Enrichment result##################
		if(!is.null(enrichedSig)){
			tabsContent <- paste(tabsContent, enrichResultTab(enrichMethod, geneSetDag), seq='\n')
			if (!is.null(geneSetDag)) {
				dagRes <- expandDag(enrichedSig[, "geneset"], geneSetDag)
				dagEdges <- dagRes$edges
				dagNodes <- getDagNodes(enrichedSig, dagRes$allNodes, geneSetDes, enrichMethod, dagColor)
				dagJson <- toJSON(unname(c(dagEdges, dagNodes)))

			}
			# tabsContent <- paste(tabsContent, enrichResultTabCategoryViz(outputHtmlFile,organism,enrichMethod,fdrMethod,enrichedSig,dNum,geneSetDag,geneSetDes,geneSetNet,outputDirectory,timeStamp,dagColor,hostName,interestingGeneMap,enrichDatabase), sep='\n')
		}

		#template <- readLines(system.file("templates/tab.mustache", package="WebGestaltR"))
		#data <- list(hasEnrichedSig=!is.null(enrichedSig), idIsEntrezGene=interestingGeneMap$standardId=="entrezgene", tabsContent=tabsContent)
		#bodyContent <- whisker.render(template, data)
		bodyContent <- tabsContent
	}else{
		###########Organism is others. No mapping information#############
		#############summary for the analysis###################
		bodyContent <- summaryDescription(timeStamp,organism,interestGeneFile,interestGene,interestGeneType,enrichMethod,enrichDatabase,enrichDatabaseFile,enrichDatabaseType,enrichDatabaseDescriptionFile,interestingGeneMap,referenceGeneList,referenceGeneFile,referenceGene,referenceGeneType,referenceSet,minNum,maxNum,sigMethod,fdrThr,topThr,fdrMethod,enrichedSig,dNum,perNum,lNum,geneSet)

		##############Enrich Result################
		if(!is.null(enrichedSig)){
			bodyContent <- paste(bodyContent, enrichResultOthers(enrichMethod,enrichedSig,geneSetDes,fdrMethod,dNum), sep='\n')
		}
		standardId <- NULL
	}
	if (is.null(enrichedSig)) {
		enrichedSig <- data.frame()
	}
	if (is.null(background)) {
		background <- data.frame()
	}
	version <- packageVersion("WebGestaltR")

	header <- readLines(system.file("templates/header.mustache", package="WebGestaltR"))
	footer <- readLines(system.file("templates/footer.mustache", package="WebGestaltR"))
	template <- readLines(system.file("templates/template.mustache", package="WebGestaltR"))
	data <- list(hostName=hostName, geneSetNet=geneSetNet, geneSetDag=geneSetDag, bodyContent=bodyContent,
				sigJson=toJSON(unname(rowSplit(enrichedSig))), insigJson=toJSON(unname(rowSplit(background))),
				dagJson=dagJson, hasGeneSetDag=!is.null(geneSetDag), version=version,
				geneTableJson=toJSON(geneTables), standardId=standardId, numAnnoRefUserId=numAnnoRefUserId,
				methodIsGsea=enrichMethod=="GSEA", hasGeneSetDes=!is.null(geneSetDes)
				)
	cat(whisker.render(template, data, partials=list(header=header, footer=footer)), file=outputHtmlFile)
}
