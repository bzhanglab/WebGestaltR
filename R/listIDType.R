listIdType <- function(organism="hsapiens",hostName="http://www.webgestalt.org/"){
	jsonData <- fromJSON(file=file.path(hostName,"data","idtypesummary.json"))
	idType <- jsonData[[organism]]
	idType <- sapply(idType,function(e){return(e$name)})
	return(idType)
}

listIDType <- function(...) {
	cat("WARNING: Function listIDType is deprecated and changed to listIdType!\n")
	return(listIdType(...))
}
