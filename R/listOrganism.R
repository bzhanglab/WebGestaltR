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
  cacheData <- cacheFile(hostName, c("summary", "idtype"))
  if (! cacheData$Succeed) {
    return(cacheData$ERROR)
  }
  jsonData <- cacheData$jsonData
	organisms <- names(jsonData)
	return(organisms)
}
