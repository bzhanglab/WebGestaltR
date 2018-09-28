createReport <- function(hostName, outputDirectory, organism="hsapiens", projectName, enrichMethod, geneSet, geneSetDes, geneSetDag, geneSetNet, interestingGeneMap, referenceGeneList, enrichedSig, geneTables, clusters, background, enrichDatabase="geneontology_Biological_Process", enrichDatabaseFile=NULL, enrichDatabaseType=NULL, enrichDatabaseDescriptionFile=NULL, interestGeneFile=NULL, interestGene=NULL, interestGeneType=NULL, collapseMethod="mean", referenceGeneFile=NULL, referenceGene=NULL, referenceGeneType=NULL, referenceSet=NULL, minNum=10, maxNum=500, fdrMethod="BH", sigMethod="fdr", fdrThr=0.05, topThr=10, dNum=20, perNum=1000, lNum=20, dagColor="binary"){

	outputHtmlFile <- file.path(outputDirectory, paste0("Project_", projectName), paste0("Report_", projectName, ".html"))

	# if hostname starts with "file://", it is used as WebGestaltReporter
	# all web assets are avaialble inside a parent directory called "assets"
	## TODO: FIXME
	if(length(grep("file://", hostName, fixed=TRUE))==1){
		#file.symlink("../assets", file.path(outputDirectory, paste("Project_",projectName,sep=""),"assets"))
		#hostName <- "assets"
		hostName <- "https://s3-us-west-2.amazonaws.com/webgestalt/assets"
	}

	numAnnoRefUserId <- NULL
	dagJson <- NULL
	if(organism!="others"){
		#####Summary Tab########
		tabsContent <- summaryDescription(projectName, organism, interestGeneFile, interestGene, interestGeneType, enrichMethod, enrichDatabase, enrichDatabaseFile, enrichDatabaseType, enrichDatabaseDescriptionFile, interestingGeneMap, referenceGeneList, referenceGeneFile, referenceGene, referenceGeneType, referenceSet, minNum, maxNum, sigMethod, fdrThr, topThr, fdrMethod, enrichedSig, dNum, perNum, lNum, geneSet)

		if (sigMethod == "fdr" && !is.null(enrichedSig) && dNum < nrow(enrichedSig)) {
			enrichedSig <- enrichedSig[1:dNum, ]
		}

		standardId <- interestingGeneMap$standardId
		if (enrichMethod == 'ORA') {
			interestGeneList <- unique(interestingGeneMap$mapped[[standardId]])
			numAnnoRefUserId <- length(intersect(interestGeneList, intersect(referenceGeneList, geneSet$gene)))
		}
		###########GOSlim summary#########################
		if(standardId=="entrezgene"){
			tabsContent <- paste(tabsContent, goSlimReport(projectName), sep='\n')
		}

		############Enrichment result##################
		if(!is.null(enrichedSig)){
			tabsContent <- paste(tabsContent, enrichResultTab(enrichMethod, geneSetDes, geneSetDag), seq='\n')
			if (!is.null(geneSetDag)) {
				dagRes <- expandDag(enrichedSig$geneSet, geneSetDag)
				dagEdges <- dagRes$edges
				dagNodes <- getDagNodes(enrichedSig, dagRes$allNodes, geneSetDes, enrichMethod, dagColor)
				dagJson <- toJSON(unname(c(dagEdges, dagNodes)))

			}
		}
		bodyContent <- tabsContent
	}else{
		###########Organism is others. No mapping information#############
		#############summary for the analysis###################
		if (enrichMethod == 'ORA') {
			numAnnoRefUserId <- length(intersect(interestingGeneMap, intersect(referenceGeneList, geneSet$gene)))
		}
		bodyContent <- summaryDescription(projectName, organism, interestGeneFile, interestGene, interestGeneType, enrichMethod, enrichDatabase, enrichDatabaseFile, enrichDatabaseType, enrichDatabaseDescriptionFile, interestingGeneMap, referenceGeneList, referenceGeneFile, referenceGene, referenceGeneType, referenceSet, minNum, maxNum, sigMethod, fdrThr, topThr, fdrMethod, enrichedSig, dNum, perNum, lNum, geneSet)

		if (sigMethod == "fdr" && !is.null(enrichedSig) && dNum < nrow(enrichedSig)) {
			enrichedSig <- enrichedSig[1:dNum, ]
		}

		##############Enrich Result################
		if(!is.null(enrichedSig)){
			bodyContent <- paste(bodyContent, enrichResultTab(enrichMethod, geneSetDes, geneSetDag), seq='\n')
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
	if (!is.null(clusters$ap)) {
		apClusterJson = toJSON(unname(clusters$ap$clusters))
		apRepJson = toJSON(unname(clusters$ap$representatives))
	} else {
		apClusterJson = "''"
		apRepJson = "''"
	}

	header <- readLines(system.file("templates/header.mustache", package="WebGestaltR"))
	footer <- readLines(system.file("templates/footer.mustache", package="WebGestaltR"))
	template <- readLines(system.file("templates/template.mustache", package="WebGestaltR"))
	data <- list(hostName=hostName, geneSetNet=geneSetNet, bodyContent=bodyContent,
				sigJson=toJSON(unname(rowSplit(enrichedSig))), insigJson=toJSON(unname(rowSplit(background))),
				dagJson=dagJson, hasGeneSetDag=!is.null(geneSetDag), version=version,
				apClusterJson=apClusterJson, apRepJson=apRepJson,
				geneTableJson=toJSON(geneTables), standardId=standardId, numAnnoRefUserId=numAnnoRefUserId,
				methodIsGsea=enrichMethod=="GSEA", hasGeneSetDes=!is.null(geneSetDes)
				)
	cat(whisker.render(template, data, partials=list(header=header, footer=footer)), file=outputHtmlFile)
}
