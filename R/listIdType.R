listIdType <- function(organism="hsapiens",hostName="http://www.webgestalt.org/"){
	response <- GET(file.path(hostName, "api", "summary", "idtype"))
	if (response$status_code != 200) {
		return(webRequestError(reponse))
	}
	jsonData <- content(response)
	idType <- jsonData[[organism]]
	idType <- sapply(idType,function(e){return(e$name)})
	return(idType)
}

listIDType <- function(...) {
	cat("WARNING: Function listIDType is deprecated and changed to listIdType!\n")
	return(listIdType(...))
}
