listOrganism <- function(hostName="http://www.webgestalt.org/"){
	jsonData <- fromJSON(file=file.path(hostName,"data","idtypesummary.json"))
	organisms <- names(jsonData)
	return(organisms)
}
