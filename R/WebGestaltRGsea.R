#' @importFrom dplyr select distinct left_join arrange %>%
#' @importFrom readr write_tsv
WebGestaltRGsea <- function(organism="hsapiens", enrichDatabase="geneontology_Biological_Process", enrichDatabaseFile=NULL, enrichDatabaseType=NULL, enrichDatabaseDescriptionFile=NULL,  interestGeneFile=NULL, interestGene=NULL, interestGeneType=NULL, collapseMethod="mean", minNum=10, maxNum=500, fdrMethod="BH", sigMethod="fdr", fdrThr=0.05, topThr=10, reportNum=20, setCoverNum=10, perNum=1000, isOutput=TRUE, outputDirectory=getwd(), projectName=NULL, dagColor="binary", nThreads=1, hostName="http://www.webgestalt.org/") {
	enrichMethod <- "GSEA"
	projectDir <- file.path(outputDirectory, paste0("Project_", projectName))

	#########Web server will input "NULL" to the R package, thus, we need to change "NULL" to NULL########
	enrichDatabase <- testNull(enrichDatabase)
	enrichDatabaseFile <- testNull(enrichDatabaseFile)
	enrichDatabaseType <- testNull(enrichDatabaseType)
	enrichDatabaseDescriptionFile <- testNull(enrichDatabaseDescriptionFile)
	interestGeneFile <- testNull(interestGeneFile)
	interestGene <- testNull(interestGene)
	interestGeneType <- testNull(interestGeneType)

	################Check parameter################
	errorTest <- parameterErrorMessage(enrichMethod=enrichMethod, organism=organism, collapseMethod=collapseMethod, minNum=minNum, maxNum=maxNum, fdrMethod=fdrMethod, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, reportNum=reportNum, perNum=perNum, isOutput=isOutput, outputDirectory=outputDirectory, dagColor=dagColor, hostName=hostName)

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
	interestingGeneMap <- loadInterestGene(organism=organism, dataType="rnk", inputGeneFile=interestGeneFile, inputGene=interestGene, geneType=interestGeneType, collapseMethod=collapseMethod, hostName=hostName, geneSet=geneSet)

	if(.hasError(interestingGeneMap)){
		return(interestingGeneMap)
	}

	if(organism=="others"){
		interestGeneList <- unique(interestingGeneMap)
	}else{
		interestStandardId <- interestingGeneMap$standardId
		interestGeneList <- interestingGeneMap$mapped %>% select(interestStandardId, .data$score) %>% distinct()
	}

	##########Create project folder##############
	if(isOutput==TRUE){
		dir.create(projectDir)

	######Summarize gene annotation based on the GOSlim###########
		if(organism!="others"){
			if(databaseStandardId=="entrezgene"){
				cat("Summarize the uploaded ID list by GO Slim data...\n")
				goSlimOutput <- file.path(projectDir, paste0("goslim_summary_", projectName))
				re <- goSlimSummary(organism=organism, geneList=interestGeneList[[interestStandardId]], outputFile=goSlimOutput, outputType="png", isOutput=isOutput, hostName=hostName)
				if(.hasError(re)){
					return(re)
				}
			}
			write_tsv(interestingGeneMap$mapped, file.path(projectDir, paste0("interestingID_mappingTable_", projectName, ".txt")))
			write(interestingGeneMap$unmapped, file.path(projectDir, paste0("interestingID_unmappedList_", projectName, ".txt")))
		}else{
			write_tsv(interestGeneList, file.path(projectDir, paste0("interestList_", projectName, ".txt")), col_names=FALSE)
		}
	}

	#############Run enrichment analysis###################
	cat("Perform the enrichment analysis...\n")

	gseaRes <- gseaEnrichment(hostName, outputDirectory, projectName, interestGeneList,
		geneSet, geneSetDes=geneSetDes, minNum=minNum, maxNum=maxNum, sigMethod=sigMethod, fdrThr=fdrThr,
		topThr=topThr, perNum=perNum, nThreads=nThreads, isOutput=isOutput
	)

	if(.hasError(gseaRes)){
		return(gseaRes)
	}

	enrichedSig <- gseaRes$enriched
	insig <- gseaRes$background


	clusters <- list()
	geneTables <- list()
	if(!is.null(enrichedSig)){
		if(!is.null(geneSetDes)){ #######Add extra description information###########
			enrichedSig <- enrichedSig %>%
				left_join(geneSetDes, by="geneSet") %>%
				select(.data$geneSet, .data$description, .data$ES, .data$NES, .data$pValue, .data$FDR, .data$link, .data$size, .data$plotPath, .data$leadingEdgeNum, .data$leadingEdgeId) %>%
				arrange(.data$FDR, .data$pValue, desc(.data$NES))
		} else {
			enrichedSig <- enrichedSig %>%
				select(.data$geneSet, .data$ES, .data$NES, .data$pValue, .data$FDR, .data$link, .data$size, .data$plotPath, .data$leadingEdgeNum, .data$leadingEdgeId) %>%
				arrange(.data$FDR, .data$pValue, desc(.data$NES))
		}

		geneTables <- getGeneTables(organism, enrichedSig, "leadingEdgeId", interestingGeneMap)
		if (organism != "others") {
			enrichedSig$link <- mapply(function(link, geneList) linkModification("GSEA", enrichDatabase, link, geneList, interestingGeneMap),
				enrichedSig$link,
				enrichedSig$leadingEdgeId
			)
		}

		if (organism != "others" && interestGeneType != interestStandardId) {
			outputEnrichedSig <- mapUserId(enrichedSig, "leadingEdgeId", interestingGeneMap)
		} else {
			outputEnrichedSig <- enrichedSig
		}

		if(isOutput==TRUE){
			write_tsv(outputEnrichedSig, file.path(projectDir, paste0("enrichment_results_", projectName, ".txt")))
			idsInSet <- sapply(enrichedSig$leadingEdgeId, strsplit, split=";")
			names(idsInSet) <- enrichedSig$geneSet
			pValue <- enrichedSig$pValue
			pValue[pValue == 0] <- .Machine$double.eps
			signedLogP <- -log(pValue) * sign(enrichedSig$ES)
			apRes <- affinityPropagation(idsInSet, signedLogP)
			wscRes <- weightedSetCover(idsInSet, 1 / signedLogP, setCoverNum, nThreads)
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
		createReport(hostName=hostName, outputDirectory=outputDirectory, organism=organism, projectName=projectName, enrichMethod=enrichMethod, geneSet=geneSet, geneSetDes=geneSetDes, geneSetDag=geneSetDag, geneSetNet=geneSetNet, interestingGeneMap=interestingGeneMap, enrichedSig=enrichedSig, background=insig, geneTables=geneTables, clusters=clusters, enrichDatabase=enrichDatabase, enrichDatabaseFile=enrichDatabaseFile, enrichDatabaseType=enrichDatabaseType, enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile, interestGeneFile=interestGeneFile, interestGene=interestGene, interestGeneType=interestGeneType, collapseMethod=collapseMethod, minNum=minNum, maxNum=maxNum, fdrMethod=fdrMethod, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, reportNum=reportNum, perNum=perNum, dagColor=dagColor)

		cwd <- getwd()
		setwd(projectDir)
		zip(paste0("Project_", projectName, ".zip"), ".", flags="-rq")
		setwd(cwd)

		cat("Results can be found in the ", projectDir, "!\n", sep="")
	}
	return(outputEnrichedSig)
}
