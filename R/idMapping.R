#' ID Mapping
#'
#' ID mapping utility with WebGestalt server.
#'
#' @inheritParams WebGestaltR
#' @param dataType Type of data, either \code{list}, \code{rnk} or \code{gmt}.
#'   Could be \code{list}, \code{rnk} or \code{matrix} for \code{idToSymbol}.
#' @param inputGeneFile The data file to be mapped.
#' @param inputGene Or the input could be given as an R object.
#'   GMT file should be read with \code{readGmt}.
#' @param sourceIdType The ID type of the data.
#' @param targetIdType The ID type of the mapped data.
#' @param mappingOutput Boolean if the mapping output is written to file.
#' @param outputFileName The output file name.
#'
#' @return A list of \code{mapped} and \code{unmapped} IDs.
#' @export
#' @aliases IDMapping
#'
idMapping <- function(organism="hsapiens", dataType="list", inputGeneFile=NULL, inputGene=NULL, sourceIdType, targetIdType=NULL, collapseMethod="mean", mappingOutput=FALSE, outputFileName="", cache=NULL, hostName="https://www.webgestalt.org/") {
	#############Check general parameters########
	errorTest <- parameterErrorMessage(organism=organism, dataType=dataType, collapseMethod=collapseMethod, hostName=hostName, mappingOutput=mappingOutput, cache=cache)

	if(!is.null(errorTest)){
		stop(errorTest)
	}

	############Check source id type#########
	errorTest <- idTypeError(idType=sourceIdType, organism=organism, hostName=hostName, cache=cache)
	if(!is.null(errorTest)){
		stop(errorTest)
	}

	##########Identify the standardId for the input ID type###########
	standardSource <- identifyStandardId(hostName=hostName, idType=sourceIdType, organism=organism, type="interest", cache=cache)

	############Check target id type#########
	if(!is.null(targetIdType)){
		errorTest <- targetIdTypeError(idType=targetIdType, organism=organism, hostName=hostName, cache=cache)
		if(!is.null(errorTest)){
			stop(errorTest)
		}else{
			standardTarget <- identifyStandardId(hostName=hostName, idType=targetIdType, organism=organism, type="interest", cache=cache)
			errorTest <- stardardDiffError(standardSource=standardSource,standardTarget=standardTarget)
			if(!is.null(errorTest)){
				stop(errorTest)
			}
		}
	}else{
		targetIdType <- standardSource
	}

	##########gene level ID Mapping##########
	if(standardSource=="entrezgene"){
		idMap <- idMappingGene(organism=organism, dataType=dataType, inputGeneFile=inputGeneFile, inputGene=inputGene, sourceIdType=sourceIdType, targetIdType=targetIdType, collapseMethod=collapseMethod, mappingOutput=mappingOutput, outputFileName=outputFileName, hostName=hostName)
	} else {
		idMap <- idMappingPhosphosite(organism=organism, dataType=dataType, inputGeneFile=inputGeneFile, inputGene=inputGene, sourceIdType=sourceIdType, targetIdType=targetIdType, collapseMethod=collapseMethod, mappingOutput=mappingOutput, outputFileName=outputFileName, hostName=hostName)
	}


	idMap$standardId <- standardSource
	return(idMap)
}

#' @export
IDMapping <- function(...) {
	warning("Function IDMapping is deprecated and changed to idMapping!\n")
	return(idMapping(...))
}
