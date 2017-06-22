listOrganism <- function(hostName="http://www.webgestalt.org/"){
	json_data <- fromJSON(file=file.path(hostName,"data","idtypesummary.json"))
	organisms <- names(json_data)
  return(organisms)
}