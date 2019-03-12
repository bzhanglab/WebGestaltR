#' @importFrom readr write_tsv
#' @importFrom dplyr left_join select arrange %>% desc
WebGestaltROra <- function(organism="hsapiens", enrichDatabase="geneontology_Biological_Process", enrichDatabaseFile=NULL, enrichDatabaseType=NULL, enrichDatabaseDescriptionFile=NULL,  interestGeneFile=NULL, interestGene=NULL, interestGeneType=NULL, collapseMethod="mean", referenceGeneFile=NULL, referenceGene=NULL, referenceGeneType=NULL, referenceSet=NULL, minNum=10, maxNum=500, fdrMethod="BH", sigMethod="fdr", fdrThr=0.05, topThr=10, reportNum=20, setCoverNum=10, isOutput=TRUE, outputDirectory=getwd(), projectName=NULL, dagColor="binary", nThreads=1, hostName="http://www.webgestalt.org/"){
	enrichMethod <- "ORA"
	projectDir <- file.path(outputDirectory, paste0("Project_", projectName))

	#########Web server will input "NULL" to the R package, thus, we need to change "NULL" to NULL########
	enrichDatabase <- testNull(enrichDatabase)
	enrichDatabaseFile <- testNull(enrichDatabaseFile)
	enrichDatabaseType <- testNull(enrichDatabaseType)
	enrichDatabaseDescriptionFile <- testNull(enrichDatabaseDescriptionFile)
	interestGeneFile <- testNull(interestGeneFile)
	interestGene <- testNull(interestGene)
	interestGeneType <- testNull(interestGeneType)
	referenceGeneFile <- testNull(referenceGeneFile)
	referenceGene <- testNull(referenceGene)
	referenceGeneType <- testNull(referenceGeneType)
	referenceSet <- testNull(referenceSet)

	################Check parameter################
	errorTest <- parameterErrorMessage(enrichMethod=enrichMethod, organism=organism, collapseMethod=collapseMethod, minNum=minNum, maxNum=maxNum, fdrMethod=fdrMethod, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, reportNum=reportNum, isOutput=isOutput, outputDirectory=outputDirectory, dagColor=dagColor, hostName=hostName)

	if(!is.null(errorTest)){
		return(errorTest)
	}

	#############Check enriched database#############
	cat("Uploading the functional categories...\n")
	enrichD <- loadGeneSet(organism=organism, enrichDatabase=enrichDatabase, enrichDatabaseFile=enrichDatabaseFile, enrichDatabaseType=enrichDatabaseType, enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile, collapseMethod=collapseMethod, hostName=hostName)
	if(.hasError(enrichD)){
		return(enrichD)
	}

	geneSet <- enrichD$geneSet
	geneSetDes <- enrichD$geneSetDes
	geneSetDag <- enrichD$geneSetDag
	geneSetNet <- enrichD$geneSetNet
	databaseStandardId <- enrichD$standardId

	###########Check input interesting gene list###############
	cat("Uploading the ID list...\n")
	interestingGeneMap <- loadInterestGene(organism=organism, dataType="list", inputGeneFile=interestGeneFile, inputGene=interestGene, geneType=interestGeneType, collapseMethod=collapseMethod, hostName=hostName, geneSet=geneSet)

	if(.hasError(interestingGeneMap)){
		return(interestingGeneMap)
	}

	if(organism=="others"){
		interestGeneList <- unique(interestingGeneMap)
	}else{
		interestStandardId <- interestingGeneMap$standardId
		interestGeneList <- unique(interestingGeneMap$mapped[[interestStandardId]])
	}

	###################load reference gene set for SEA method##############
	cat("Uploading the reference list...\n")
	referenceGeneList <- loadReferenceGene(organism=organism, referenceGeneFile=referenceGeneFile, referenceGene=referenceGene, referenceGeneType=referenceGeneType, referenceSet=referenceSet, collapseMethod=collapseMethod, hostName=hostName, geneSet=geneSet, interestGeneList=interestGeneList)

	if(.hasError(referenceGeneList)){
		return(referenceGeneList)
	}

	##########Create project folder##############
	if(isOutput==TRUE){
		dir.create(projectDir)

	######Summarize gene annotation based on the GOSlim###########
		if(organism!="others"){
			if(databaseStandardId=="entrezgene"){
				cat("Summarize the uploaded ID list by GO Slim data...\n")
				goSlimOutput <- file.path(projectDir, paste0("goslim_summary_", projectName))
				re <- goSlimSummary(organism=organism, geneList=interestGeneList, outputFile=goSlimOutput, outputType="png", isOutput=isOutput, hostName=hostName)
				if(.hasError(re)){
					return(re)
				}
			}
			write_tsv(interestingGeneMap$mapped, file.path(projectDir, paste0("interestingID_mappingTable_", projectName , ".txt")))
			write(interestingGeneMap$unmapped, file.path(projectDir, paste0("interestingID_unmappedList_", projectName, ".txt")))
		}else{
			write(interestGeneList, file.path(projectDir, paste0("interestList_", projectName, ".txt")))
		}
	}

	#############Run enrichment analysis###################
	cat("Perform the enrichment analysis...\n")

	oraRes <- oraEnrichment(interestGeneList,referenceGeneList,geneSet,minNum=minNum,maxNum=maxNum,fdrMethod=fdrMethod,sigMethod=sigMethod,fdrThr=fdrThr,topThr=topThr)

	if(.hasError(oraRes)){
		return(oraRes)
	}

	enrichedSig <- oraRes$enriched
	insig <- oraRes$background

	clusters <- list()
	geneTables <- list()

	if(!is.null(enrichedSig)){
		if(!is.null(geneSetDes)){ #######Add extra description information###########
			enrichedSig <- enrichedSig %>%
				left_join(geneSetDes, by="geneSet") %>%
				select(.data$geneSet, .data$description, .data$link, .data$C, .data$O, .data$E, .data$R, .data$pValue, .data$FDR, .data$overlapId) %>%
				arrange(.data$FDR, .data$pValue, desc(.data$C))
		} else {
			enrichedSig <- enrichedSig %>%
				select(.data$geneSet, .data$link, .data$C, .data$O, .data$E, .data$R, .data$pValue, .data$FDR, .data$overlapId) %>%
				arrange(.data$FDR, .data$pValue, desc(.data$C))
		}

		geneTables <- getGeneTables(organism, enrichedSig, "overlapId", interestingGeneMap)
		if (organism != "others") {
			enrichedSig$link <- mapply(function(link, geneList) linkModification("ORA", enrichDatabase, link, geneList, interestingGeneMap),
				enrichedSig$link,
				enrichedSig$overlapId
			)
		}

		if (organism != "others" && interestGeneType != interestStandardId) {
			outputEnrichedSig <- mapUserId(enrichedSig, "overlapId", interestingGeneMap)
		} else {
			outputEnrichedSig <- enrichedSig
		}

		if(isOutput==TRUE){
			write_tsv(outputEnrichedSig, file.path(projectDir, paste0("enrichment_results_", projectName, ".txt")))
			idsInSet <- sapply(enrichedSig$overlapId, strsplit, split=";")
			names(idsInSet) <- enrichedSig$geneSet
			minusLogP <- -log(enrichedSig$pValue)
			minusLogP[minusLogP == Inf] <- -log(.Machine$double.eps)
			apRes <- affinityPropagation(idsInSet, minusLogP)
			wscRes <- weightedSetCover(idsInSet, 1 / minusLogP, setCoverNum, nThreads)
			if (!is.null(apRes)) {
				writeLines(sapply(apRes$clusters, paste, collapse="\t"), file.path(projectDir, paste0("enriched_geneset_ap_clusters_", projectName, ".txt")))
			} else {
				apRes <- NULL
			}
			clusters$ap <- apRes
			if (!is.null(wscRes$topSets)) {
				writeLines(c(paste0("# Coverage: ", wscRes$coverage), wscRes$topSets), file.path(projectDir, paste0("enriched_geneset_wsc_topsets_", projectName, ".txt")))
				clusters$wsc <- list(representatives=wscRes$topSets,  coverage=wscRes$coverage)
			} else {
				clusters$wsc <- NULL
			}
		}
	}

	if(isOutput==TRUE){

	##############Create report##################
		cat("Generate the final report...\n")
		createReport(hostName=hostName, outputDirectory=outputDirectory, organism=organism, projectName=projectName, enrichMethod=enrichMethod, geneSet=geneSet, geneSetDes=geneSetDes, geneSetDag=geneSetDag, geneSetNet=geneSetNet, interestingGeneMap=interestingGeneMap, referenceGeneList=referenceGeneList, enrichedSig=enrichedSig, background=insig, geneTables=geneTables, clusters=clusters, enrichDatabase=enrichDatabase, enrichDatabaseFile=enrichDatabaseFile, enrichDatabaseType=enrichDatabaseType, enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile, interestGeneFile=interestGeneFile, interestGene=interestGene, interestGeneType=interestGeneType, collapseMethod=collapseMethod, referenceGeneFile=referenceGeneFile, referenceGene=referenceGene, referenceGeneType=referenceGeneType, referenceSet=referenceSet, minNum=minNum, maxNum=maxNum, fdrMethod=fdrMethod, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, reportNum=reportNum, dagColor=dagColor)

		cwd <- getwd()
		setwd(projectDir)
		zip(paste0("Project_", projectName, ".zip"), ".", flags="-rq")
		setwd(cwd)

		cat("Results can be found in the ", projectDir, "!\n", sep="")
	}

	return(outputEnrichedSig)
}
