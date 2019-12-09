#' @importFrom dplyr select distinct left_join arrange %>% mutate
#' @importFrom readr write_tsv
WebGestaltRGsea <- function(organism="hsapiens", enrichDatabase=NULL, enrichDatabaseFile=NULL, enrichDatabaseType=NULL, enrichDatabaseDescriptionFile=NULL,  interestGeneFile=NULL, interestGene=NULL, interestGeneType=NULL, collapseMethod="mean", minNum=10, maxNum=500, fdrMethod="BH", sigMethod="fdr", fdrThr=0.05, topThr=10, reportNum=20, setCoverNum=10, perNum=1000, p=1, isOutput=TRUE, outputDirectory=getwd(), projectName=NULL, dagColor="binary", saveRawGseaResult=FALSE, plotFormat="png", nThreads=1, cache=NULL, hostName="http://www.webgestalt.org/") {
	enrichMethod <- "GSEA"
	projectDir <- file.path(outputDirectory, paste0("Project_", projectName))

	######### Web server will input "NULL" to the R package, thus, we need to change "NULL" to NULL ########
	enrichDatabase <- testNull(enrichDatabase)
	enrichDatabaseFile <- testNull(enrichDatabaseFile)
	enrichDatabaseType <- testNull(enrichDatabaseType)
	enrichDatabaseDescriptionFile <- testNull(enrichDatabaseDescriptionFile)
	interestGeneFile <- testNull(interestGeneFile)
	interestGene <- testNull(interestGene)
	interestGeneType <- testNull(interestGeneType)

	################ Check parameter ################
	errorTest <- parameterErrorMessage(enrichMethod=enrichMethod, organism=organism, collapseMethod=collapseMethod, minNum=minNum, maxNum=maxNum, fdrMethod=fdrMethod, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, reportNum=reportNum, perNum=perNum, isOutput=isOutput, outputDirectory=outputDirectory, dagColor=dagColor, hostName=hostName, cache=cache)

	if(!is.null(errorTest)){
		stop(errorTest)
	}

	############# Check enriched database #############
	cat("Loading the functional categories...\n")
	enrichD <- loadGeneSet(organism=organism, enrichDatabase=enrichDatabase, enrichDatabaseFile=enrichDatabaseFile, enrichDatabaseType=enrichDatabaseType, enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile, cache=cache, hostName=hostName)

	geneSet <- enrichD$geneSet
	geneSetDes <- enrichD$geneSetDes
	geneSetDag <- enrichD$geneSetDag
	geneSetNet <- enrichD$geneSetNet
	databaseStandardId <- enrichD$standardId
	rm(enrichD)

	########### Check input interesting gene list ###############
	cat("Loading the ID list...\n")
	interestingGeneMap <- loadInterestGene(organism=organism, dataType="rnk", inputGeneFile=interestGeneFile, inputGene=interestGene, geneType=interestGeneType, collapseMethod=collapseMethod, cache=cache, hostName=hostName, geneSet=geneSet)

	if (organism == "others") {
		interestGeneList <- unique(interestingGeneMap)
	} else {
		interestStandardId <- interestingGeneMap$standardId
		interestGeneList <- interestingGeneMap$mapped %>% select(interestStandardId, .data$score) %>% distinct()
	}

	########## Create project folder ##############
	if (isOutput) {
		dir.create(projectDir)

	###### Summarize gene annotation based on the GOSlim ###########
		if (organism != "others") {
			if (databaseStandardId == "entrezgene") {
				cat("Summarizing the uploaded ID list by GO Slim data...\n")
				goSlimOutput <- file.path(projectDir, paste0("goslim_summary_", projectName))
				re <- goSlimSummary(organism=organism, geneList=interestGeneList[[interestStandardId]], outputFile=goSlimOutput, outputType="png", isOutput=isOutput, cache=cache, hostName=hostName)
			}
			write_tsv(interestingGeneMap$mapped, file.path(projectDir, paste0("interestingID_mappingTable_", projectName, ".txt")))
			write(interestingGeneMap$unmapped, file.path(projectDir, paste0("interestingID_unmappedList_", projectName, ".txt")))
		} else {
			write_tsv(interestGeneList, file.path(projectDir, paste0("interestList_", projectName, ".txt")), col_names=FALSE)
		}
	}

	############# Run enrichment analysis ###################
	cat("Performing the enrichment analysis...\n")

	gseaRes <- gseaEnrichment(hostName, outputDirectory, projectName, interestGeneList,
		geneSet, geneSetDes=geneSetDes, minNum=minNum, maxNum=maxNum, sigMethod=sigMethod, fdrThr=fdrThr,
		topThr=topThr, perNum=perNum, p=p, nThreads=nThreads, saveRawGseaResult=saveRawGseaResult, plotFormat=plotFormat, isOutput=isOutput
	)
	if (is.null(gseaRes)) {
		return(NULL)
	}

	enrichedSig <- gseaRes$enriched
	insig <- gseaRes$background


	clusters <- list()
	geneTables <- list()
	if (!is.null(enrichedSig)) {
		if (!is.null(geneSetDes)) { ####### Add extra description information ###########
			enrichedSig <- enrichedSig %>%
				left_join(geneSetDes, by="geneSet") %>%
				select(.data$geneSet, .data$description, .data$link, .data$enrichmentScore, .data$normalizedEnrichmentScore, .data$pValue, .data$FDR, .data$size, .data$plotPath, .data$leadingEdgeNum, .data$leadingEdgeId) %>%
				arrange(.data$FDR, .data$pValue, desc(.data$normalizedEnrichmentScore)) %>%
				mutate(description=ifelse(is.na(.data$description), "", .data$description))
		} else {
			enrichedSig <- enrichedSig %>%
				select(.data$geneSet, .data$link, .data$enrichmentScore, .data$normalizedEnrichmentScore, .data$pValue, .data$FDR, .data$size, .data$plotPath, .data$leadingEdgeNum, .data$leadingEdgeId) %>%
				arrange(.data$FDR, .data$pValue, desc(.data$normalizedEnrichmentScore))
		}

		geneTables <- getGeneTables(organism, enrichedSig, "leadingEdgeId", interestingGeneMap)
		if (organism != "others") {
			enrichedSig$link <- mapply(function(link, geneList) linkModification("GSEA", link, geneList, interestingGeneMap),
				enrichedSig$link,
				enrichedSig$leadingEdgeId
			)
		}

		if ("database" %in% colnames(geneSet)) {
			# add source database for multiple databases
			enrichedSig <- enrichedSig %>% left_join(unique(geneSet[, c("geneSet", "database")]), by="geneSet")
		}

		if (organism != "others" && interestGeneType != interestStandardId) {
			outputEnrichedSig <- mapUserId(enrichedSig, "leadingEdgeId", interestingGeneMap)
		} else {
			outputEnrichedSig <- enrichedSig
		}

		if (isOutput) {
			write_tsv(outputEnrichedSig, file.path(projectDir, paste0("enrichment_results_", projectName, ".txt")))
			idsInSet <- sapply(enrichedSig$leadingEdgeId, strsplit, split=";")
			names(idsInSet) <- enrichedSig$geneSet
			pValue <- enrichedSig$pValue
			pValue[pValue == 0] <- .Machine$double.eps
			signedLogP <- -log(pValue) * sign(enrichedSig$enrichmentScore)
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

	if (isOutput) {
		############## Create report ##################
		cat("Generate the final report...\n")
		createReport(hostName=hostName, outputDirectory=outputDirectory, organism=organism, projectName=projectName, enrichMethod=enrichMethod, geneSet=geneSet, geneSetDes=geneSetDes, geneSetDag=geneSetDag, geneSetNet=geneSetNet, interestingGeneMap=interestingGeneMap, enrichedSig=enrichedSig, background=insig, geneTables=geneTables, clusters=clusters, enrichDatabase=enrichDatabase, enrichDatabaseFile=enrichDatabaseFile, enrichDatabaseType=enrichDatabaseType, enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile, interestGeneFile=interestGeneFile, interestGene=interestGene, interestGeneType=interestGeneType, collapseMethod=collapseMethod, minNum=minNum, maxNum=maxNum, fdrMethod=fdrMethod, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, reportNum=reportNum, perNum=perNum, p=p, dagColor=dagColor)

		cwd <- getwd()
		setwd(projectDir)
		zip(paste0("Project_", projectName, ".zip"), ".", flags="-rq")
		setwd(cwd)

		cat("Results can be found in the ", projectDir, "!\n", sep="")
	}
	return(outputEnrichedSig)
}
