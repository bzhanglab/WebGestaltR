#' List Reference Sets
#'
#' List available reference sets for the given organism on WebGestalt server.
#'
#' @inheritParams WebGestaltR
#'
#' @return A list of reference sets.
#'
#' @importFrom httr GET content
#' @export
#'
listReferenceSet <- function(organism="hsapiens",hostName="http://www.webgestalt.org/"){
	if (startsWith(hostName, "file://")) {
		jsonData <- fromJSON(file=removeFileProtocol(file.path(hostName, "referencesetsummary.json")))
	} else {
		response <- GET(file.path(hostName, "api", "summary", "referenceset"))
		if (response$status_code != 200) {
			return(webRequestError(response))
		}
		jsonData <- content(response)
	}
	idType <- jsonData[[organism]]
	idType <- sapply(idType,function(e){return(e$name)})
	return(idType)
}
