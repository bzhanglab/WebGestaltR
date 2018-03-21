listReferenceSet <- function(organism="hsapiens",hostName="http://www.webgestalt.org/"){
	jsonData <- fromJSON(file=file.path(hostName,"data","referenceSetsummary.json"))
	idType <- jsonData[[organism]]
	idType <- sapply(idType,function(e){return(e$name)})
	return(idType)
}
