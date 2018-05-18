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
		interestGeneList <- unique(interestingGeneMap[,1])
	}else{
		interestStandardId <- interestingGeneMap$standardId
		interestGeneList <- unique(interestingGeneMap$mapped[,interestStandardId])
	}

	###################load reference gene set for SEA method##############
	cat("Uploading the reference list...\n")
	referenceGeneList <- NULL

	referenceGeneList <- loadReferenceGene(organism=organism, referenceGeneFile=referenceGeneFile, referenceGene=referenceGene, referenceGeneType=referenceGeneType, referenceSet=referenceSet, collapseMethod=collapseMethod, hostName=hostName, geneSet=geneSet, interestGeneList=interestGeneList)

	if(.hasError(referenceGeneList)){
		return(referenceGeneList)
	}

	##########Create project folder##############
	if(isOutput==TRUE){
		dir.create(projectDir)

	######Summary gene annotation based on the GOSlim###########
		if(organism!="others"){
			if(databaseStandardId=="entrezgene"){
				cat("Summary the uploaded ID list by GO Slim data...\n")
				goSlimOutput <- file.path(projectDir,paste("goslim_summary_",timeStamp,sep=""))
				re <- goSlimSummary(organism=organism,geneList=interestGeneList,outputFile=goSlimOutput,outputType="png",hostName=hostName)
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

	enrichedSig <- oraEnrichment(interestGeneList,referenceGeneList,geneSet,minNum=minNum,maxNum=maxNum,fdrMethod=fdrMethod,sigMethod=sigMethod,fdrThr=fdrThr,topThr=topThr)

	if(.hasError(enrichedSig)){
		return(enrichedSig)
	}

	if(!is.null(enrichedSig)){
		if(!is.null(geneSetDes)){ #######Add extra description information###########
			colnames(geneSetDes) <- c("geneset","description")
			enrichedSig <- merge(x=enrichedSig,y=geneSetDes,by="geneset",all.x=TRUE)
			enrichedSig <- enrichedSig[,c(1,ncol(enrichedSig),2:(ncol(enrichedSig)-1))]
			enrichedSig <- enrichedSig[order(enrichedSig[,"FDR"],enrichedSig[,"PValue"]),]
		}

		if(organism!="others" && interestGeneType!=interestStandardId){
			enrichedSig <- mapUserId(enrichedSig,"overlapID",interestingGeneMap)  ###mapUserId function is in the enrichmentResultProcess_component file
		}

		if(isOutput==TRUE){
			write.table(enrichedSig,file=file.path(projectDir,paste("enrichment_results_",timeStamp,".txt",sep="")),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
		}
	}

	if(isOutput==TRUE){

	##############Create report##################
		cat("Generate the final report...\n")
		createReport(hostName=hostName,outputDirectory=outputDirectory,organism=organism,timeStamp=timeStamp,enrichMethod=enrichMethod,geneSet=geneSet,geneSetDes=geneSetDes,geneSetDag=geneSetDag,geneSetNet=geneSetNet,interestingGeneMap=interestingGeneMap,referenceGeneList=referenceGeneList,enrichedSig=enrichedSig,enrichDatabase=enrichDatabase,enrichDatabaseFile=enrichDatabaseFile,enrichDatabaseType=enrichDatabaseType,enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile,interestGeneFile=interestGeneFile,interestGene=interestGene,interestGeneType=interestGeneType,collapseMethod=collapseMethod,referenceGeneFile=referenceGeneFile,referenceGene=referenceGene,referenceGeneType=referenceGeneType,referenceSet=referenceSet,minNum=minNum,maxNum=maxNum,fdrMethod=fdrMethod,sigMethod=sigMethod,fdrThr=fdrThr,topThr=topThr,dNum=dNum,dagColor=dagColor)

		zip(file.path(projectDir, paste0(projectDir, ".zip")), projectDir, flags="-jrq")

		cat("Results can be found in the ",projectDir,"!",sep="")
	}

	return(enrichedSig)
}
