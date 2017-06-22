listReferenceSet <- function(organism="hsapiens",hostName="http://www.webgestalt.org/"){
	json_data <- fromJSON(file=file.path(hostName,"data","referenceSetsummary.json"))
	idtype <- json_data[[organism]]
	idtype <- unlist(lapply(idtype,function(e){return(e$name)}))
	return(idtype)
}