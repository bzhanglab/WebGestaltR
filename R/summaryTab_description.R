summaryTab_description <- function(timeStamp,organism,interestGeneFile,interestGene,interestGeneType,enrichMethod,enrichDatabase,enrichDatabaseFile,enrichDatabaseType,enrichDatabaseDescriptionFile,interestingGeneMap,referenceGene_List,referenceGeneFile,referenceGene,referenceGeneType,referenceSet,minNum,maxNum,sigMethod,fdrThr,topThr,fdrMethod,enrichedSig,dNum,perNum,lNum,geneSet){
	if(enrichMethod=="ORA"){
		methodSpecificContent <- specificParameterSummaryOra(organism,referenceGene_List,geneSet,referenceGeneFile,referenceGene,referenceGeneType,referenceSet,minNum,maxNum,sigMethod,fdrThr,topThr,fdrMethod,enrichedSig,dNum,interestingGeneMap)
	}

	if(enrichMethod=="GSEA"){
		methodSpecificContent <- specificParameterSummaryGsea(organism,interestingGeneMap,geneSet,minNum,maxNum,sigMethod,fdrThr,topThr,perNum,lNum,enrichedSig,dNum)
	}

	standardId <- interestingGeneMap$standardId
	template <- readLines(system.file("inst/templates/summaryTab.mustache", package="WebGestaltR"))
	data <- list(timeStamp=timeStamp, enrichMethod=enrichMethod, organism=organism, organismIsOthers=organism=="others",
		enrichDatabase=enrichDatabase, enrichDatabaseIsOthers=enrichDatabase=="others", enrichDatabaseFile=enrichDatabaseFile,
		enrichDatabaseType=enrichDatabaseType, enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile,
		hasEnrichDatabaseDescriptioFile=!is.null(enrichDatabaseDescriptionFile), hasInterestGeneFile=!is.null(interestGeneFile),
		interestGeneFileBase=basename(interestGeneFile), interestGeneType=interestGeneType,
		numUserId=nrow(interestingGeneMap$mapped)+length(interestingGeneMap$unmapped),
		numMappedUserId=nrow(interestingGeneMap$mapped), numUniqueMappedId=length(unique(interestingGeneMap$mapped[,standardId])),
		numUnmappedUserId=length(interestingGeneMap$unmapped), idIsEntrezGene=standardId[["staIDs"]]=="entrezgene", standardId=standardId,
		methodSpecificContent=methodSpecificContent
		)

	return(whisker.render(template, data))
}


specificParameterSummaryOra <- function(organism,referenceGene_List,geneSet,referenceGeneFile,referenceGene,referenceGeneType,referenceSet,minNum,maxNum,sigMethod,fdrThr,topThr,fdrMethod,enrichedSig,dNum,interestingGeneMap){
	organismIsOthers <- organism == "others"
	if(!organismIsOthers){
		standardId <- interestingGeneMap$standardId
		interestGene_List <- unique(interestingGeneMap$mapped[,standardId])
		numAnnoRefUserId <- length(intersect(interestGene_List,intersect(referenceGene_List,geneSet[,3])))
	}else{  ###for others
		standardId <- NULL
		interestGene_List <- unique(interestingGeneMap[,1])
		numAnnoRefUserId <- NULL

	}
	numAnnoRefId <- length(intersect(referenceGene_List,geneSet[,3]))
	hasEnrichedSig <- !is.null(enrichedSig)
	if (hasEnrichedSig) {
		numEnrichedSig <- nrow(enrichedSig)
		showAll <- dNum>=numEnrichedSig
	} else {
		numEnrichedSig <- NULL
		showAll <- NULL
	}
	data <- list(organismIsOthers=organismIsOthers, numUniqueUserId=length(interestGene_List), standardId=standardId,
		numAnnoRefUserId=numAnnoRefUserId, hasRefGeneFile=!is.null(referenceGeneFile), referenceGeneFile=referenceGeneFile,
		referenceGeneType=referenceGeneType, hasRefGene=!is.null(referenceGene), referenceSet=referenceSet,
		numRefGene=length(referenceGene_List), numAnnoRefId=length(intersect(referenceGene_List,geneSet[,3])), minNum=minNum,
		maxNum=maxNum, fdrMethod=fdrMethod, methodIsFdr=sigMethod=="fdr", methodIsTop=sigMethod=="top", fdrThr=fdrThr,
		topThr=topThr, hasEnrichedSig=hasEnrichedSig, showAll=showAll, numEnrichedSig=numEnrichedSig
		)
	template <- readLines(system.file("inst/templates/summaryTabOra.mustache", package="WebGestaltR"))
	return(whisker.render(template, data))
}


specificParameterSummaryGsea <- function(organism,interestingGeneMap,geneSet,minNum,maxNum,sigMethod,fdrThr,topThr,perNum,lNum,enrichedSig,dNum){
	organismIsOthers <- organism == "others"
	if(!organismIsOthers){
		standardId <- interestingGeneMap$standardId
		interestGene_List <- unique(interestingGeneMap$mapped[,standardId])
		numUniqueUserId <- nrow(interestGene_List)
		numAnnoUserId <- length(intersect(interestGene_List,geneSet[,3]))
	} else {
		standardId <- NULL
		interestGene_List <- unique(interestingGeneMap[,1])
		numUniqueUserId <- length(interestGene_List)
		numAnnoUserId <- NULL
	}

	hasEnrichedSig <- !is.null(enrichedSig)
	data <- list(organismIsOthers=organismIsOthers, numUniqueUserId=numUniqueUserId, standardId=standardId, numAnnoUserId=numAnnoUserId,
		minNum=minNum, maxNum=maxNum, methodIsFdr=sigMethod=="fdr", methodIsTop=sigMethod=="top", fdrThr=fdrThr, topThr=topThr,
		perNum=perNum, lNum=lNum, hasEnrichedSig=hasEnrichedSig
		)

	if(hasEnrichedSig){
		x <- enrichedSig[enrichedSig[,"NES"]>0,]
		y <- enrichedSig[enrichedSig[,"NES"]<0,]
		data$numPosRel <- nrow(x)
		data$numNegRel <- nrow(y)
		data$isPosRel <- data$numPosRel>0
		data$isNegRel <- data$numNegRel>0
		data$showAllPos <- dNum>=data$numPosRel
		data$showAllNeg <- dNum>=data$numNegRel
	}

	template <- readLines(system.file("inst/templates/summaryTabGsea.mustache", package="WebGestaltR"))
	return(whisker.render(template, data))
}
