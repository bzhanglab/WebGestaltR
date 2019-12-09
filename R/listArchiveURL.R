#' List WebGestalt Servers
#'
#' List available WebGestalt servers.
#'
#'
#' @return A data frame of available servers.
#'
#' @importFrom readr read_tsv
#' @export
#' @aliases listArchiveURL
#'
listArchiveUrl <- function(){
	archiveUrl <- read_tsv("http://www.webgestalt.org/archiveURL.txt", col_names=FALSE)
	return(archiveUrl)
}

#' @export
listArchiveURL <- function(...) {
	warning("Function listArchiveURL is deprecated and changed to listArchiveUrl!\n")
	return(listArchiveUrl(...))
}
