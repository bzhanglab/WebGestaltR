listGeneSet <- function(organism="hsapiens",hostName="http://www.webgestalt.org/"){

	json_data <- fromJSON(file=file.path(hostName,"data","genesetsummary.json"))
	ids <- json_data[[organism]]
	name1 <- names(ids)
	idList <- data.frame(name="",description="",idtype="",stringsAsFactors=F)
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
