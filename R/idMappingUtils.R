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
identifyStandardId <- function(hostName,idType,organism,type){
	if(type=="interest"){
		cacheData <- cacheFile(hostName, c("summary", "idtype"))
	}
	if(type=="reference"){
	  cacheData <- cacheFile(hostName, c("summary", "referenceset"))
	}
	if (! cacheData$Succeed) {
		return(cacheData$ERROR)
	}
	jsonData <- cacheData$jsonData
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
