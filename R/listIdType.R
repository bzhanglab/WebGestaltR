#' List ID Types
#'
#' List supported ID types for the given organism on WebGestalt server.
#'
#' @inheritParams WebGestaltR
#'
#' @return A list of supported gene sets.
#'
#' @importFrom httr content
#' @export
#' @aliases listIDType
#'
listIdType <- function(organism="hsapiens", hostName="https://www.webgestalt.org/", cache=NULL) {
	if (startsWith(hostName, "file://")) {
        jsonData <- fromJSON(removeFileProtocol(file.path(hostName, "idtypesummary.json")))
		idType <- jsonData[[organism]]
		idType <- idType$name
	} else {
		response <- cacheUrl(file.path(hostName, "api", "summary", "idtype"), cache)
		if (response$status_code != 200) {
			return(webRequestError(response))
		}
		jsonData <- content(response)
        idType <- jsonData[[organism]]
        idType <- sapply(idType,function(e){return(e$name)})
    }
	return(idType)
}

#' @export
listIDType <- function(...) {
	warning("Function listIDType is deprecated and changed to listIdType!\n")
	return(listIdType(...))
}
