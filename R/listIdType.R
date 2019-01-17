#' List ID Types
#'
#' List supported ID types for the given organism on WebGestalt server.
#'
#' @inheritParams WebGestaltR
#'
#' @return A list of supported gene sets.
#'
#' @importFrom httr GET content
#' @export
#' @aliases listIDType
#'
listIdType <- function(organism="hsapiens",hostName="http://www.webgestalt.org/"){
	response <- GET(file.path(hostName, "api", "summary", "idtype"))
	if (response$status_code != 200) {
		return(webRequestError(response))
	}
	jsonData <- content(response)
	idType <- jsonData[[organism]]
	idType <- sapply(idType,function(e){return(e$name)})
	return(idType)
}

#' @export
listIDType <- function(...) {
	cat("WARNING: Function listIDType is deprecated and changed to listIdType!\n")
	return(listIdType(...))
}
