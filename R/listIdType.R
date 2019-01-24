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
  cacheData <- cacheFile(hostName, c("summary", "idtype"))
  if (! cacheData$Succeed) {
    return(cacheData$ERROR)
  }
  jsonData <- cacheData$jsonData
	idType <- jsonData[[organism]]
	idType <- sapply(idType,function(e){return(e$name)})
	return(idType)
}

#' @export
listIDType <- function(...) {
	cat("WARNING: Function listIDType is deprecated and changed to listIdType!\n")
	return(listIdType(...))
}
