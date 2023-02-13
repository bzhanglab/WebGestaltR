#' List Reference Sets
#'
#' List available reference sets for the given organism on WebGestalt server.
#'
#' @inheritParams WebGestaltR
#'
#' @return A list of reference sets.
#'
#' @importFrom httr content
#' @export
#'
listReferenceSet <- function(organism="hsapiens", hostName="https://www.webgestalt.org/", cache=NULL) {
	if (startsWith(hostName, "file://")) {
		jsonData <- fromJSON(removeFileProtocol(file.path(hostName, "referencesetsummary.json")))
		idType <- jsonData[[organism]]
		idType <- idType$name
	} else {
		response <- cacheUrl(file.path(hostName, "api", "summary", "referenceset"), cache)
		if (response$status_code != 200) {
			return(webRequestError(response))
		}
		jsonData <- content(response)
		idType <- jsonData[[organism]]
		idType <- sapply(idType,function(e){return(e$name)})
	}
	return(idType)
}
