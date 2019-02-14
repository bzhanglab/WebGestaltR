#' List Organisms
#'
#' List supported organisms on WebGestalt server.
#'
#' @inheritParams WebGestaltR
#'
#' @return A list of supported organisms.
#'
#' @importFrom httr GET content
#' @importFrom rjson fromJSON
#' @export
#'
listOrganism <- function(hostName="http://www.webgestalt.org/"){
	if (startsWith(hostName, "file://")) {
		jsonData <- fromJSON(file=removeFileProtocol(file.path(hostName, "idtypesummary.json")))
	} else {
		response <- GET(file.path(hostName, "api", "summary", "idtype"))
		if (response$status_code != 200) {
			return(webRequestError(response))
		}
		jsonData <- content(response)
	}
	organisms <- names(jsonData)
	return(organisms)
}
