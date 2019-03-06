loadInterestGene <- function(organism="hsapiens", dataType="list", inputGeneFile=NULL, inputGene=NULL, geneType="entrezgene", collapseMethod="mean", hostName="http://www.webgestalt.org/", geneSet){
	if(is.null(inputGeneFile) && is.null(inputGene)){
		return(interestGeneError(type="empty"))
	}else{
		if(organism!="others"){
			if(is.null(geneType)){
				return(interestGeneError(type="emptyType"))
			}else{
				mapRe <- .uploadGeneExistingOrganism(organism=organism, dataType=dataType, inputGeneFile=inputGeneFile, inputGene=inputGene, geneType=geneType, collapseMethod=collapseMethod, geneSet=geneSet, hostName=hostName)
				if(.hasError(mapRe)){
					return(mapRe)
				}
			}
		}else{
			mapRe <- .uploadGeneOthers(dataType=dataType,inputGeneFile=inputGeneFile,inputGene=inputGene,geneSet=geneSet)
			if(.hasError(mapRe)){
				return(mapRe)
			}
		}
	}

	#if organism is not others, the function will return a mapping result with mapped and unmapped list
	#if organism is others, the function will return a matrix with gene list
	return(mapRe)
}

#' @importFrom httr GET content
#' @importFrom readr read_tsv
loadReferenceGene <- function(organism="hsapiens", referenceGeneFile=NULL, referenceGene=NULL, referenceGeneType="entrezgene", referenceSet=NULL, collapseMethod="mean", hostName="http://www.webgestalt.org/", geneSet, interestGeneList) {
	referenceGeneList <- NULL
	referenceGeneMap <- NULL

	if(is.null(referenceGeneFile) && is.null(referenceGene) && is.null(referenceSet)){
		return(referenceGeneError(type="empty"))
	}else{
		if(organism!="others"){
			if(!is.null(referenceGeneFile) || !is.null(referenceGene)){
				if(is.null(referenceGeneType)){
					return(referenceGeneError(type="emptyType"))
				}else{
					mapRe <- .uploadGeneExistingOrganism(organism=organism, dataType="list", inputGeneFile=referenceGeneFile, inputGene=referenceGene, geneType=referenceGeneType, collapseMethod=collapseMethod, geneSet=geneSet, hostName=hostName)
					if(.hasError(mapRe)){
						return(mapRe)
					}
					geneStandardId <- identifyStandardId(hostName=hostName,idType=referenceGeneType,organism=organism,type="interest")
					referenceGeneList <- mapRe$mapped[[geneStandardId]]
				}
			}else{ ###referenceGeneFile and referenceGene are both NULL. But referenceSet is not NULL
				refS <- listReferenceSet(organism=organism,hostName=hostName)
				if(length(which(refS==referenceSet))==0){
					return(referenceGeneError(type="existingRef"))
				}
				refStandardId <- identifyStandardId(hostName=hostName,idType=referenceSet,organism=organism,type="reference")
				if (startsWith(hostName, "file://")) {
					# Getting data from local directory in the old way
					refPath <- removeFileProtocol(file.path(hostName, "reference", paste0(paste(organism, referenceSet, refStandardId, sep="_"), ".table")))
					referenceGeneList <- read_tsv(refPath, col_names=FALSE, col_types="c-")[[1]]
				} else {
					response <- GET(file.path(hostName, "api", "reference"), query=list(organism=organism, referenceSet=referenceSet, standardId=refStandardId))
					if (response$status_code != 200) {
						return(webRequestError(response))
					}
					# API now just returns one single column
					referenceGeneList <- read_tsv(content(response), col_names=FALSE, col_types="c")[[1]]
				}
			}
		}else{ ##For other organisms
			if(!is.null(referenceGeneFile) || !is.null(referenceGene)){
				referenceGeneList <- .uploadGeneOthers(dataType="list",inputGeneFile=referenceGeneFile,inputGene=referenceGene,geneSet=geneSet)
				if(.hasError(referenceGeneList)){
					return(referenceGeneList)
				}
				referenceGeneList <- unique(referenceGeneList)
			}else{
				return(referenceGeneError(type="empty"))
			}
		}
	}

	##compare interest gene list and reference gene list
	if(length(intersect(interestGeneList, intersect(referenceGeneList, geneSet$gene)))==0){
		return(referenceGeneError(type="interestEmpty"))
	}
	return(referenceGeneList)
}


#' @importFrom dplyr filter
.uploadGeneExistingOrganism <- function(organism, dataType, inputGeneFile, inputGene, geneType, collapseMethod, geneSet, hostName) {
	geneMap <- idMapping(organism=organism, dataType=dataType, inputGeneFile=inputGeneFile, inputGene=inputGene, sourceIdType=geneType, targetIdType=NULL, collapseMethod=collapseMethod, mappingOutput=FALSE, hostName=hostName)

	if(.hasError(geneMap)){
		return(geneMap)
	}

	#gene_standardId <- identifyStandardId(hostName=hostName,idtype=geneType,organism=organism,type="interest")  ##identifyStandardId in idMappingComponent.R
	#if(gene_standardId!=databaseStandardId){  ###the standardId of the input genes should be the same with the standardarId of the functional database
#	return(interestGeneError(type="unmatch"))
#}

	geneMapMappedList <- geneMap$mapped
	standardId <- geneMap$standardId

	geneList <- as.character(unique(geneMapMappedList[[standardId]]))
	ov <- intersect(geneList, geneSet$gene)

	if(length(ov)==0){
		return(interestGeneError(type="unannotated"))
	}

	###Because if all genes are annotated to only one category, GSEA will return the error, we need to avoid this error by reporting the error in the R#
	geneSets <- unique((filter(geneSet, .data$gene %in% geneList))[["geneSet"]])
	if (length(geneSets) == 1) {
		return(interestGeneError(type="onlyOne"))
	}
	return(geneMap)
}


#' @importFrom dplyr filter
.uploadGeneOthers <- function(dataType,inputGeneFile,inputGene,geneSet){
	inputGene <- formatCheck(dataType=dataType,inputGeneFile=inputGeneFile,inputGene=inputGene)
	if(.hasError(inputGene)){
		return(inputGene)
	}

	if (dataType == "list") {
		geneList = inputGene
	} else if (dataType == "rnk") {
		geneList = inputGene$gene
	}

	ov <- intersect(geneList, geneSet$gene)
	if(length(ov)==0){
		return(interestGeneError(type="unannotated"))
	}

	###Because if all genes are annotated to only one category, GSEA will return the error, we need to avoid this error by reporting the error in the R#
	geneSets <- unique((filter(geneSet, .data$gene %in% geneList))[["geneSet"]])
	if (length(geneSets) == 1) {
		return(interestGeneError(type="onlyOne"))
	}
	return(inputGene)
}
