idMappingInput <- function(dataType="list",inputGeneFile,inputGene){
	if(dataType=="gmt"){
		if(!is.null(inputGeneFile)){
			inputGene <- readGmt(inputGeneFile)
			return(inputGene)
		}else{
			stop(gmtFormatError("empty"))
		}
	}else{
		inputGene <- formatCheck(dataType=dataType,inputGeneFile=inputGeneFile,inputGene=inputGene)
		return(inputGene)
	}
}

#' @importFrom httr content
#' @importFrom jsonlite fromJSON
identifyStandardId <- function(hostName, idType, organism, type, cache) {
	if (startsWith(hostName, "file://")) {
		if (type=="interest") {
			summaryPath <- removeFileProtocol(file.path(hostName, "idtypesummary.json"))
		}
		if (type=="reference") {
			summaryPath <- removeFileProtocol(file.path(hostName, "referencesetsummary.json"))
		}
		jsonData <- fromJSON(summaryPath)
		idTypes <- jsonData[[organism]]
		names <- idTypes$name
		standardIds <- idTypes$type
	} else {
		if (type=="interest") {
			response <- cacheUrl(file.path(hostName, "api", "summary", "idtype"), cache)
		}
		if (type=="reference") {
			response <- cacheUrl(file.path(hostName, "api", "summary", "referenceset"), cache)
		}
		if (response$status_code != 200) {
			stop(webRequestError(response))
		}
		jsonData <- content(response)
		idTypes <- jsonData[[organism]]
		names <- unlist(lapply(idTypes, function(e) return(e$name)))
		standardIds <- unlist(lapply(idTypes,function(e) return(e$type)))
	}
	idTypes <- data.frame(name=names, standardId=standardIds, stringsAsFactors=FALSE)
	return(filter(idTypes, .data$name == idType)[[1, "standardId"]])
}

#' @importFrom dplyr select distinct %>%
#' @importFrom readr write_tsv
idMappingOutput <- function(outputFileName, mappingList, unmappedList, dataType, sourceIdType, targetIdType) {
	if (length(unmappedList)>0) {
		write(unmappedList, paste0(outputFileName, "_unmappedList.txt"))
	}
	if (dataType == "list" | dataType == "rnk") {
		dataType = "txt"
	}
	fileName <- paste0(outputFileName, "_mappedList_from_", sourceIdType, "_to_", targetIdType, ".", dataType)
	if(dataType=="gmt"){
		genes <- tapply(mappingList[[targetIdType]], mappingList$geneSet, paste, collapse="\t")
		gmtDf <- mappingList %>% select(.data$geneSet, .data$link) %>% distinct()
		gmtDf$genes = genes[gmtDf$geneSet]
		write_tsv(gmtDf, fileName, col_names=FALSE)
	}else{
		write_tsv(mappingList, fileName)
	}
}

replace_null <- function(x) {
	lapply(x, function(x) {
		if (is.list(x)) {
			replace_null(x)
		} else {
		if (is.null(x)) NA else(x)
		}
	})
}
