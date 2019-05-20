urlToFile <- function(dataUrl) {
	result <- sub("^http://", "", dataUrl)
	result <- sub("^https://", "", result)
	result <- gsub("\\?[^?]+?=", "_", result)
	result <- gsub("&[^&]+?=", "_", result)
	result <- gsub("[:/.]", "_", result)
	result <- gsub("_+", "_", result, fixed=TRUE)
	return(result)
}

#' cacheUrl
#'
#' Get data from a URL or cache and optionally save in cache for reuse
#'
#' @param dataUrl The URL of data
#' @param cache The cache directory. Defaults to \code{NULL} and fetches data from server.
#' @param query The list of queries passed on to httr methods
#'
#' @return response object from httr request
#'
#' @importFrom httr GET
#' @keywords internal
#'
cacheUrl <- function(dataUrl, cache=NULL, query=NULL) {
	if (!is.null(cache)) {
		dir.create(cache, showWarnings=FALSE)
		if (!is.null(query)) {
			localFilePrefix <- urlToFile(paste0(dataUrl, "_", paste0(query, collapse="_")))
		} else {
			localFilePrefix <- urlToFile(dataUrl)
		}
		localFile <- file.path(cache, paste0(localFilePrefix, ".rds"))
	}
	if (!is.null(cache) && file.exists(localFile)) {
		#cat("Reading from cache: ", localFile, "\n")
		response <- readRDS(localFile)
	} else {
		#cat("Reading from server: ", dataUrl, "\n")
		if (!is.null(query)) {
			response <- GET(dataUrl, query=query)
		} else {
			response <- GET(dataUrl)
		}
		if (response$status_code != 200) {
			return(response)
		}
		if (!is.null(cache)) {
			saveRDS(response, localFile)
		}
	}
	return(response)
}
