idMappingInput <- function(dataType="list",inputGeneFile,inputGene){
	if(dataType=="gmt"){
		if(!is.null(inputGeneFile)){
			inputGene <- readGmt(inputGeneFile)
			return(inputGene)
		}else{
			return(gmtFormatError("empty"))
		}
	}else{
		inputGene <- formatCheck(dataType=dataType,inputGeneFile=inputGeneFile,inputGene=inputGene)
		return(inputGene)
	}
}

#' @importFrom httr GET content
#' @importFrom rjson fromJSON
identifyStandardId <- function(hostName,idType,organism,type){
	if (startsWith(hostName, "file://")) {
		if (type=="interest") {
			summaryPath <- removeFileProtocol(file.path(hostName, "idtypesummary.json"))
		}
		if (type=="reference") {
			summaryPath <- removeFileProtocol(file.path(hostName, "referencesetsummary.json"))
		}
		jsonData <- fromJSON(file=summaryPath)
	} else {
		if (type=="interest") {
			response <- GET(file.path(hostName, "api", "summary", "idtype"))
		}
		if (type=="reference") {
			response <- GET(file.path(hostName, "api", "summary", "referenceset"))
		}
		if (response$status_code != 200) {
			return(webRequestError(response))
		}
		jsonData <- content(response)
	}
	idTypes <- jsonData[[organism]]
	names <- unlist(lapply(idTypes, function(e) return(e$name)))
	standardIds <- unlist(lapply(idTypes,function(e) return(e$type)))
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
