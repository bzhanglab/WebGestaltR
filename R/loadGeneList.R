loadInterestGene <- function(organism="hsapiens", dataType="list", inputGeneFile=NULL, inputGene=NULL, geneType="entrezgene", collapseMethod="mean", cache=NULL, hostName="https://www.webgestalt.org/", geneSet) {
	if (is.null(inputGeneFile) && is.null(inputGene)) {
		stop(interestGeneError(type="empty"))
	} else {
		if (organism!="others") {
			if (is.null(geneType)) {
				stop(interestGeneError(type="emptyType"))
			} else {
				mapRe <- .uploadGeneExistingOrganism(organism=organism, dataType=dataType, inputGeneFile=inputGeneFile, inputGene=inputGene, geneType=geneType, collapseMethod=collapseMethod, geneSet=geneSet, cache=cache, hostName=hostName)
			}
		} else {
			mapRe <- .uploadGeneOthers(dataType=dataType, inputGeneFile=inputGeneFile, inputGene=inputGene, geneSet=geneSet)
		}
	}

	#if organism is not others, the function will return a mapping result with mapped and unmapped list
	#if organism is others, the function will return a matrix with gene list
	return(mapRe)
}

#' @importFrom httr content
#' @importFrom readr read_tsv
loadReferenceGene <- function(organism="hsapiens", referenceGeneFile=NULL, referenceGene=NULL, referenceGeneType="entrezgene", referenceSet=NULL, collapseMethod="mean", hostName="https://www.webgestalt.org/", geneSet, interestGeneList, cache=NULL) {
	referenceGeneList <- NULL
	referenceGeneMap <- NULL

	if (is.null(referenceGeneFile) && is.null(referenceGene) && is.null(referenceSet)) {
		stop(referenceGeneError(type="empty"))
	} else {
		if (organism!="others") {
			if (!is.null(referenceGeneFile) || !is.null(referenceGene)) {
				if (is.null(referenceGeneType)) {
					stop(referenceGeneError(type="emptyType"))
				} else {
					mapRe <- .uploadGeneExistingOrganism(organism=organism, dataType="list", inputGeneFile=referenceGeneFile, inputGene=referenceGene, geneType=referenceGeneType, collapseMethod=collapseMethod, geneSet=geneSet, cache=cache, hostName=hostName)
					geneStandardId <- identifyStandardId(hostName=hostName, idType=referenceGeneType, organism=organism, type="interest", cache=cache)
					referenceGeneList <- mapRe$mapped[[geneStandardId]]
				}
			} else { ### referenceGeneFile and referenceGene are both NULL. But referenceSet is not NULL
				refS <- listReferenceSet(organism=organism, hostName=hostName, cache=cache)
				if (length(which(refS==referenceSet))==0) {
					stop(referenceGeneError(type="existingRef"))
				}
				refStandardId <- identifyStandardId(hostName=hostName, idType=referenceSet, organism=organism, type="reference", cache=cache)
				if (startsWith(hostName, "file://")) {
					# Getting data from local directory in the old way
					refPath <- removeFileProtocol(file.path(hostName, "reference", paste0(paste(organism, referenceSet, refStandardId, sep="_"), ".table")))
					referenceGeneList <- read_tsv(refPath, col_names=FALSE, col_types="c-")[[1]]
				} else {
					response <- cacheUrl(file.path(hostName, "api", "reference"), cache=cache, query=list(organism=organism, referenceSet=referenceSet, standardId=refStandardId))
					if (response$status_code != 200) {
						stop(webRequestError(response))
					}
					# API now just returns one single column
					referenceGeneList <- read_tsv(content(response), col_names=FALSE, col_types="c")[[1]]
				}
			}
		} else { ## For other organisms
			if (!is.null(referenceGeneFile) || !is.null(referenceGene)) {
				referenceGeneList <- .uploadGeneOthers(dataType="list", inputGeneFile=referenceGeneFile, inputGene=referenceGene, geneSet=geneSet)
				referenceGeneList <- unique(referenceGeneList)
			} else {
				stop(referenceGeneError(type="empty"))
			}
		}
	}

	## compare interest gene list and reference gene list
	if (length(intersect(interestGeneList, intersect(referenceGeneList, geneSet$gene)))==0) {
		stop(referenceGeneError(type="interestEmpty"))
	}
	return(referenceGeneList)
}


#' @importFrom dplyr filter
.uploadGeneExistingOrganism <- function(organism, dataType, inputGeneFile, inputGene, geneType, collapseMethod, geneSet, cache, hostName) {
	geneMap <- idMapping(organism=organism, dataType=dataType, inputGeneFile=inputGeneFile, inputGene=inputGene, sourceIdType=geneType, targetIdType=NULL, collapseMethod=collapseMethod, mappingOutput=FALSE, cache=cache, hostName=hostName)

	#gene_standardId <- identifyStandardId(hostName=hostName,idtype=geneType,organism=organism,type="interest")  ##identifyStandardId in idMappingComponent.R
	#if(gene_standardId!=databaseStandardId){  ###the standardId of the input genes should be the same with the standardarId of the functional database
#	return(interestGeneError(type="unmatch"))
#}

	geneMapMappedList <- geneMap$mapped
	standardId <- geneMap$standardId

	geneList <- as.character(unique(geneMapMappedList[[standardId]]))
	ov <- intersect(geneList, geneSet$gene)

	if (length(ov)==0) {
		stop(interestGeneError(type="unannotated"))
	}

	###Because if all genes are annotated to only one category, GSEA will return the error, we need to avoid this error by reporting the error in the R#
	geneSets <- unique((filter(geneSet, .data$gene %in% geneList))[["geneSet"]])
	if (length(geneSets) == 1) {
		stop(interestGeneError(type="onlyOne"))
	}
	return(geneMap)
}


#' @importFrom dplyr filter
.uploadGeneOthers <- function(dataType,inputGeneFile,inputGene,geneSet){
	inputGene <- formatCheck(dataType=dataType,inputGeneFile=inputGeneFile,inputGene=inputGene)

	if (dataType == "list") {
		geneList = inputGene
	} else if (dataType == "rnk") {
		geneList = inputGene$gene
	}

	ov <- intersect(geneList, geneSet$gene)
	if (length(ov)==0) {
		stop(interestGeneError(type="unannotated"))
	}

	###Because if all genes are annotated to only one category, GSEA will return the error, we need to avoid this error by reporting the error in the R#
	geneSets <- unique((filter(geneSet, .data$gene %in% geneList))[["geneSet"]])
	if (length(geneSets) == 1) {
		stop(interestGeneError(type="onlyOne"))
	}
	return(inputGene)
}
