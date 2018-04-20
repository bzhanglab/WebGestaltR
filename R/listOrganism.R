listOrganism <- function(hostName="http://www.webgestalt.org/"){
	response <- GET(file.path(hostName, "api", "summary", "idtype"))
	if (response$status_code != 200) {
		return(webRequestError(reponse))
	}
	jsonData <- content(response)
	organisms <- names(jsonData)
	return(organisms)
}
