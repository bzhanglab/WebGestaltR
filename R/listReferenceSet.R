listReferenceSet <- function(organism="hsapiens",hostName="http://www.webgestalt.org/"){
	response <- GET(file.path(hostName, "api", "summary", "referenceset"))
	if (response$status_code != 200) {
		return(webRequestError(response))
	}
	jsonData <- content(response)
	idType <- jsonData[[organism]]
	idType <- sapply(idType,function(e){return(e$name)})
	return(idType)
}
