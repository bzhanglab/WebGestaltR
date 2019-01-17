#' List Organisms
#'
#' List supported organisms on WebGestalt server.
#'
#' @inheritParams WebGestaltR
#'
#' @return A list of supported organisms.
#'
#' @importFrom httr GET content
#' @export
#'
listOrganism <- function(hostName="http://www.webgestalt.org/"){
	response <- GET(file.path(hostName, "api", "summary", "idtype"))
	if (response$status_code != 200) {
		return(webRequestError(response))
	}
	jsonData <- content(response)
	organisms <- names(jsonData)
	return(organisms)
}
