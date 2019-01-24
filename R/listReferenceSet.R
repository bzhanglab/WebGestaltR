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
  cacheData <- cacheFile(hostName, c("summary", "referenceset"))
  if (! cacheData$Succeed) {
    return(cacheData$ERROR)
  }
  jsonData <- cacheData$jsonData
	idType <- jsonData[[organism]]
	idType <- sapply(idType,function(e){return(e$name)})
	return(idType)
}
