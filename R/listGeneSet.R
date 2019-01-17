#' List Gene Sets
#'
#' List available gene sets for the given organism on WebGestalt server.
#'
#' @inheritParams WebGestaltR
#'
#' @return A data frame of available gene sets.
#'
#' @importFrom httr GET content
#' @export
#'
listGeneSet <- function(organism="hsapiens",hostName="http://www.webgestalt.org/"){
	response <- GET(file.path(hostName, "api", "summary", "geneset"))
	if (response$status_code != 200) {
		return(webRequestError(response))
	}
	jsonData <- content(response)
	ids <- jsonData[[organism]]
	name1 <- names(ids)
	idList <- data.frame(name="",description="",idType="",stringsAsFactors=F)
	di <- 1
	for(i in c(1:length(name1))){
		x <- .getName(ids[[i]])
		if(!(is.null(x) || (is.list(x) && length(x)==0))){
			y <- .getDescription(ids[[i]])
			z <- .getIdType(ids[[i]])
			x <- paste(name1[i],"_",x,sep="")
			idList[di:(di+length(x)-1),1] <- x
			idList[di:(di+length(x)-1),2] <- y
			idList[di:(di+length(x)-1),3] <- z
			di <- di+length(x)
		}
	}
	return(idList)
}

.getName <- function(ids){
	return(sapply(ids,function(e){return(e$name)}))
}

.getDescription <- function(ids){
	return(sapply(ids,function(e){return(e$description)}))
}

.getIdType <- function(ids){
	return(sapply(ids,function(e){return(e$idtype)}))
}
