WebGestaltRGsea <- function(organism="hsapiens", enrichDatabase="geneontology_Biological_Process", enrichDatabaseFile=NULL, enrichDatabaseType=NULL, enrichDatabaseDescriptionFile=NULL,  interestGeneFile=NULL, interestGene=NULL, interestGeneType=NULL, collapseMethod="mean", minNum=10, maxNum=500, fdrMethod="BH", sigMethod="fdr", fdrThr=0.05, topThr=10, dNum=20, perNum=1000, lNum=20, isOutput=TRUE, outputDirectory=getwd(), projectName=NULL, dagColor="binary", hostName="http://www.webgestalt.org/") {
	enrichMethod <- "GSEA"

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

	################Check parameter################
	errorTest <- parameterErrorMessage(enrichMethod=enrichMethod, organism=organism, collapseMethod=collapseMethod, minNum=minNum, maxNum=maxNum, fdrMethod=fdrMethod, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, dNum=dNum, perNum=perNum, lNum=lNum, isOutput=isOutput, outputDirectory=outputDirectory, dagColor=dagColor, hostName=hostName)

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
		interestGeneList <- interestingGeneMap$mapped %>% select(interestStandardId, score) %>% distinct()
	}

	##########Create project folder##############
	if(isOutput==TRUE){
		dir.create(projectDir)

	######Summary gene annotation based on the GOSlim###########
		if(organism!="others"){
			if(databaseStandardId=="entrezgene"){
				cat("Summary the uploaded ID list by GO Slim data...\n")
				goSlimOutput <- file.path(projectDir,paste("goslim_summary_",timeStamp,sep=""))
				re <- goSlimSummary(organism=organism,geneList=interestGeneList[[interestStandardId]],outputFile=goSlimOutput,outputType="png",hostName=hostName)
				if(.hasError(re)){
					return(re)
				}
			}
			write.table(interestingGeneMap$mapped,file=file.path(projectDir,paste("interestingID_Mappingtable_",timeStamp,".txt",sep="")),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
			write.table(interestingGeneMap$unmapped,file=file.path(projectDir,paste("interestingID_unmappedList_",timeStamp,".txt",sep="")),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
		}else{
			write.table(interestGeneList,file=file.path(projectDir,paste("interestList_",timeStamp,".txt",sep="")),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
		}
	}

	#############Run enrichment analysis###################
	cat("Perform the enrichment analysis...\n")

	gseaRes <- gseaEnrichment(hostName, outputDirectory, timeStamp, interestGeneList,
		geneSet, minNum=minNum, maxNum=maxNum, sigMethod=sigMethod, fdrThr=fdrThr,
		topThr=topThr, perNum=perNum, lNum=lNum, isOutput=isOutput
	)

	if(.hasError(gseaRes)){
		return(gseaRes)
	}

	enrichedSig <- gseaRes$enriched
	insig <- gseaRes$background

	geneTables <- list()
	if(!is.null(enrichedSig)){
		if(!is.null(geneSetDes)){ #######Add extra description information###########
			colnames(geneSetDes) <- c("geneset","description")
			enrichedSig <- enrichedSig %>%
				left_join(geneSetDes, by="geneset") %>%
				select(geneset, description, ES, NES, PValue, FDR, link, Size, plotPath, leadingEdgeNum, leadingEdgeID) %>%
				arrange(FDR, PValue)
		}

		if (organism != "others") {
			geneTables <- getGeneTables(enrichedSig, "leadingEdgeID", interestingGeneMap)
			enrichedSig$link <- mapply(function(link, geneList) linkModification(enrichDatabase, link, geneList, interestingGeneMap),
				enrichedSig$link,
				enrichedSig$leadingEdgeID
			)
		}

		if(isOutput==TRUE){
			if(organism!="others" && interestGeneType!=interestStandardId){
				outputEnrichedSig <- mapUserId(enrichedSig,"leadingEdgeID",interestingGeneMap)  ###mapUserId function is in the enrichmentResultProcess_component file
			}
			write.table(outputEnrichedSig,file=file.path(projectDir,paste("enrichment_results_",timeStamp,".txt",sep="")),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
		}
	}

	if(isOutput==TRUE){
		##############Create report##################
		cat("Generate the final report...\n")
		createReport(hostName=hostName, outputDirectory=outputDirectory, organism=organism, timeStamp=timeStamp, enrichMethod=enrichMethod, geneSet=geneSet, geneSetDes=geneSetDes, geneSetDag=geneSetDag, geneSetNet=geneSetNet, interestingGeneMap=interestingGeneMap, enrichedSig=enrichedSig, geneTables=geneTables, background=insig, enrichDatabase=enrichDatabase, enrichDatabaseFile=enrichDatabaseFile, enrichDatabaseType=enrichDatabaseType, enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile, interestGeneFile=interestGeneFile, interestGene=interestGene, interestGeneType=interestGeneType, collapseMethod=collapseMethod, minNum=minNum, maxNum=maxNum, fdrMethod=fdrMethod, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, dNum=dNum, perNum=perNum, lNum=lNum, dagColor=dagColor)

		cwd <- getwd()
		setwd(projectDir)
		zip(paste0("Project_", timeStamp, ".zip"), ".", flags="-rq")
		setwd(cwd)

		cat("Results can be found in the ",projectDir,"!",sep="")
	}
	return(enrichedSig)
}
