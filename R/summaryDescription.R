#' summaryDescription
#'
#' Render job summary section
#'
#' @importFrom whisker whisker.render
#' @keywords internal
#'
summaryDescription <- function(projectName, organism, interestGeneFile, interestGene, interestGeneType, enrichMethod, enrichDatabase, enrichDatabaseFile, enrichDatabaseType, enrichDatabaseDescriptionFile, interestingGeneMap, referenceGeneList, referenceGeneFile, referenceGene, referenceGeneType, referenceSet, minNum, maxNum, sigMethod, fdrThr, topThr, fdrMethod, enrichedSig, reportNum, perNum, geneSet, repAdded, hostName) {
	if(enrichMethod=="ORA"){
		methodSpecificContent <- specificParameterSummaryOra(organism, referenceGeneList, geneSet, referenceGeneFile, referenceGene, referenceGeneType, referenceSet, minNum, maxNum, sigMethod, fdrThr, topThr, fdrMethod, enrichedSig, reportNum, repAdded, interestingGeneMap)
	}

	if(enrichMethod=="GSEA"){
		methodSpecificContent <- specificParameterSummaryGsea(organism, interestingGeneMap, geneSet, minNum, maxNum, sigMethod, fdrThr, topThr, perNum, enrichedSig, reportNum, repAdded)
	}

	template <- readLines(system.file("templates/summary.mustache", package="WebGestaltR"))
	if (organism != "others") {
		standardId <- unname(interestingGeneMap$standardId)
		data <- list(projectName=projectName, enrichMethod=enrichMethod, organism=organism, organismIsOthers=FALSE,
			enrichDatabase=enrichDatabase, enrichDatabaseIsOthers=enrichDatabase=="others", enrichDatabaseFile=enrichDatabaseFile,
			enrichDatabaseType=enrichDatabaseType, enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile,
			hasEnrichDatabaseDescriptioFile=!is.null(enrichDatabaseDescriptionFile), hasInterestGeneFile=!is.null(interestGeneFile),
			interestGeneFileBase=ifelse(is.null(interestGeneFile), "", basename(interestGeneFile)), interestGeneType=interestGeneType,
			numUserId=nrow(interestingGeneMap$mapped)+length(interestingGeneMap$unmapped),
			numMappedUserId=nrow(interestingGeneMap$mapped), numUniqueMappedId=length(unique(interestingGeneMap$mapped[[standardId]])),
			numUnmappedUserId=length(interestingGeneMap$unmapped), idIsEntrezGene=standardId=="entrezgene", standardId=standardId,
			methodSpecificContent=methodSpecificContent, hostName=hostName
		)
	}else {
		data <- list(projectName=projectName, enrichMethod=enrichMethod, organism=organism, organismIsOthers=TRUE,
			enrichDatabase=enrichDatabase, enrichDatabaseIsOthers=enrichDatabase=="others", enrichDatabaseFile=enrichDatabaseFile,
			enrichDatabaseType=enrichDatabaseType, enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile,
			hasEnrichDatabaseDescriptioFile=!is.null(enrichDatabaseDescriptionFile), hasInterestGeneFile=!is.null(interestGeneFile),
			interestGeneFileBase=ifelse(is.null(interestGeneFile), "", basename(interestGeneFile)), interestGeneType=interestGeneType,
			idIsEntrezGene=FALSE, methodSpecificContent=methodSpecificContent
		)
	}

	return(whisker.render(template, data))
}


#' specificParameterSummaryOra
#'
#' Render job summary section of ORA specific parameters
#'
#' @importFrom whisker whisker.render
#'
#' @keywords internal
#'
specificParameterSummaryOra <- function(organism, referenceGeneList, geneSet, referenceGeneFile, referenceGene, referenceGeneType, referenceSet, minNum, maxNum, sigMethod, fdrThr, topThr, fdrMethod, enrichedSig, reportNum, repAdded, interestingGeneMap) {
	organismIsOthers <- organism == "others"
	if(!organismIsOthers){
		standardId <- interestingGeneMap$standardId
		interestGeneList <- unique(interestingGeneMap$mapped[[standardId]])
		numAnnoRefUserId <- length(intersect(interestGeneList, intersect(referenceGeneList, geneSet$gene)))
	}else{  ###for others
		standardId <- NULL
		interestGeneList <- unique(interestingGeneMap)
		numAnnoRefUserId <- NULL

	}
	numAnnoRefId <- length(intersect(referenceGeneList, geneSet$gene))
	hasEnrichedSig <- !is.null(enrichedSig)
	if (hasEnrichedSig) {
		numEnrichedSig <- nrow(enrichedSig)
		showAll <- reportNum>=numEnrichedSig
	} else {
		numEnrichedSig <- NULL
		showAll <- NULL
	}
	data <- list(organismIsOthers=organismIsOthers, numUniqueUserId=length(interestGeneList), standardId=standardId,
		numAnnoRefUserId=numAnnoRefUserId, hasRefGeneFile=!is.null(referenceGeneFile), referenceGeneFile=referenceGeneFile,
		referenceGeneType=referenceGeneType, hasRefGene=!is.null(referenceGene), referenceSet=referenceSet,
		numRefGene=length(referenceGeneList), numAnnoRefId=length(intersect(referenceGeneList, geneSet$gene)), minNum=minNum,
		maxNum=maxNum, fdrMethod=fdrMethod, methodIsFdr=sigMethod=="fdr", methodIsTop=sigMethod=="top", fdrThr=fdrThr,
		topThr=topThr, hasEnrichedSig=hasEnrichedSig, showAll=showAll, reportNum=reportNum, hasRepAdded=repAdded, numEnrichedSig=numEnrichedSig
		)
	template <- readLines(system.file("templates/summaryOra.mustache", package="WebGestaltR"))
	return(whisker.render(template, data))
}

#' specificParameterSummaryGsea
#'
#' Render job summary section of GSEA specific parameters
#'
#' @importFrom whisker whisker.render
#' @importFrom dplyr filter
#'
#' @keywords internal
#'
specificParameterSummaryGsea <- function(organism, interestingGeneMap, geneSet, minNum, maxNum, sigMethod, fdrThr, topThr, perNum, enrichedSig, reportNum, repAdded) {
	organismIsOthers <- organism == "others"
	if(!organismIsOthers){
		standardId <- interestingGeneMap$standardId
		interestGeneList <- unique(interestingGeneMap$mapped[[standardId]])
		numUniqueUserId <- length(interestGeneList)
		numAnnoUserId <- length(intersect(interestGeneList, geneSet$gene))
	} else {
		standardId <- NULL
		interestGeneList <- unique(interestingGeneMap)
		numUniqueUserId <- nrow(interestGeneList)
		numAnnoUserId <- NULL
	}

	hasEnrichedSig <- !is.null(enrichedSig)
	data <- list(organismIsOthers=organismIsOthers, numUniqueUserId=numUniqueUserId, standardId=standardId, numAnnoUserId=numAnnoUserId,
		minNum=minNum, maxNum=maxNum, methodIsFdr=sigMethod=="fdr", methodIsTop=sigMethod=="top", fdrThr=fdrThr, topThr=topThr,
		perNum=perNum, reportNum=reportNum, hasRepAdded=repAdded, hasEnrichedSig=hasEnrichedSig
		)

	if(hasEnrichedSig){
		data$numPosRel <- nrow(filter(enrichedSig, .data$NES>0))
		data$numNegRel <- nrow(filter(enrichedSig, .data$NES<0))
		data$isPosRel <- data$numPosRel>0
		data$isNegRel <- data$numNegRel>0
		data$showAll <- reportNum >= nrow(enrichedSig)
	}

	template <- readLines(system.file("templates/summaryGsea.mustache", package="WebGestaltR"))
	return(whisker.render(template, data))
}
