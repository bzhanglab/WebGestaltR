#' List Organisms
#'
#' List supported organisms on WebGestalt server.
#'
#' @inheritParams WebGestaltR
#'
#' @return A list of supported organisms.
#'
#' @importFrom httr content
#' @importFrom jsonlite fromJSON
#' @export
#'
listOrganism <- function(hostName="http://www.webgestalt.org/", cache=NULL) {
	if (startsWith(hostName, "file://")) {
		jsonData <- fromJSON(file=removeFileProtocol(file.path(hostName, "idtypesummary.json")))
	} else {
		response <- cacheUrl(file.path(hostName, "api", "summary", "idtype"), cache)
		if (response$status_code != 200) {
			return(webRequestError(response))
		}
		jsonData <- content(response)
	}
	organisms <- names(jsonData)
	return(organisms)
}
