WebGestaltROra <- function(organism="hsapiens", enrichDatabase="geneontology_Biological_Process", enrichDatabaseFile=NULL, enrichDatabaseType=NULL, enrichDatabaseDescriptionFile=NULL,  interestGeneFile=NULL, interestGene=NULL, interestGeneType=NULL, collapseMethod="mean", referenceGeneFile=NULL, referenceGene=NULL, referenceGeneType=NULL, referenceSet=NULL, minNum=10, maxNum=500, fdrMethod="BH", sigMethod="fdr", fdrThr=0.05, topThr=10, dNum=20, isOutput=TRUE, outputDirectory=getwd(), projectName=NULL, dagColor="binary", hostName="http://www.webgestalt.org/"){
	enrichMethod <- "ORA"

	if(is.null(projectName)){
		timeStamp <- as.character(as.integer(Sys.time()))
	}else{
		timeStamp <- projectName
	}

	projectDir <- file.path(outputDirectory,paste("Project_",timeStamp,sep=""))

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
	errorTest <- parameterErrorMessage(enrichMethod=enrichMethod, organism=organism, collapseMethod=collapseMethod, minNum=minNum, maxNum=maxNum, fdrMethod=fdrMethod, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, dNum=dNum, isOutput=isOutput, outputDirectory=outputDirectory, dagColor=dagColor, hostName=hostName)

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
				goSlimOutput <- file.path(projectDir,paste("goslim_summary_",timeStamp,sep=""))
				re <- goSlimSummary(organism=organism,geneList=interestGeneList,outputFile=goSlimOutput,outputType="png",hostName=hostName)
				if(.hasError(re)){
					return(re)
				}
			}
			write_tsv(interestingGeneMap$mapped, file.path(projectDir, paste0("interestingID_mappingTable_", timeStamp , ".txt")))
			write(interestingGeneMap$unmapped, file.path(projectDir, paste0("interestingID_unmappedList_", timeStamp, ".txt")))
		}else{
			write(interestGeneList, file.path(projectDir, paste0("interestList_", timeStamp, ".txt")))
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



	if(!is.null(enrichedSig)){
		if(!is.null(geneSetDes)){ #######Add extra description information###########
			colnames(geneSetDes) <- c("geneset","description")
			enrichedSig <- enrichedSig %>%
				left_join(geneSetDes, by="geneset") %>%
				select(geneset, description, link, C, O, E, R, PValue, FDR, overlapID) %>%
				arrange(FDR, PValue)
		}


		geneTables <- getGeneTables(organism, enrichedSig, "overlapID", interestingGeneMap)
		if (organism != "others") {
			enrichedSig$link <- mapply(function(link, geneList) linkModification(enrichDatabase, link, geneList, interestingGeneMap),
				enrichedSig$link,
				enrichedSig$overlapID
			)
		}

		if(isOutput==TRUE){
			if(organism!="others" && interestGeneType!=interestStandardId){
				outputEnrichedSig <- mapUserId(enrichedSig, "overlapID", interestingGeneMap)
			} else {
				outputEnrichedSig <- enrichedSig
			}
			write_tsv(outputEnrichedSig, file.path(projectDir, paste0("enrichment_results_", timeStamp, ".txt")))
			idsInSet <- sapply(enrichedSig$overlapID, strsplit, split=";")
			names(idsInSet) <- enrichedSig$geneset
			apRes <- affinityPropagation(idsInSet, enrichedSig$PValue)
			writeLines(sapply(apRes$clusters, paste, collapse="\t"), file.path(projectDir, paste0("enriched_geneset_ap_clusters_", timeStamp, ".txt")))
		}
	}

	if(isOutput==TRUE){

	##############Create report##################
		cat("Generate the final report...\n")
		createReport(hostName=hostName, outputDirectory=outputDirectory, organism=organism, timeStamp=timeStamp, enrichMethod=enrichMethod, geneSet=geneSet, geneSetDes=geneSetDes, geneSetDag=geneSetDag, geneSetNet=geneSetNet, interestingGeneMap=interestingGeneMap, referenceGeneList=referenceGeneList, enrichedSig=enrichedSig, background=insig, geneTables=geneTables, enrichDatabase=enrichDatabase, enrichDatabaseFile=enrichDatabaseFile, enrichDatabaseType=enrichDatabaseType, enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile, interestGeneFile=interestGeneFile, interestGene=interestGene, interestGeneType=interestGeneType, collapseMethod=collapseMethod, referenceGeneFile=referenceGeneFile, referenceGene=referenceGene, referenceGeneType=referenceGeneType, referenceSet=referenceSet, minNum=minNum, maxNum=maxNum, fdrMethod=fdrMethod, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, dNum=dNum, dagColor=dagColor)

		cwd <- getwd()
		setwd(projectDir)
		zip(paste0("Project_", timeStamp, ".zip"), ".", flags="-rq")
		setwd(cwd)

		cat("Results can be found in the ", projectDir, "!\n", sep="")
	}

	return(enrichedSig)
}
