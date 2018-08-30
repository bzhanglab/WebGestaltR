idMappingInput <- function(dataType="list",inputGeneFile,inputGene){
	if(dataType=="gmt"){
		if(!is.null(inputGeneFile)){
			inputGene <- readGmt(inputGeneFile)
			return(inputGene)
		}else{
			return(gmtFormatError("empty"))
		}
	}else{
		inputGene <- formatCheck(dataType=dataType,inputGeneFile=inputGeneFile,inputGene=inputGene)
		return(inputGene)
	}
}


identifyStandardId <- function(hostName,idType,organism,type){
	if(type=="interest"){
		response <- GET(file.path(hostName, "api", "summary", "idtype"))
	}
	if(type=="reference"){
		response <- GET(file.path(hostName, "api", "summary", "referenceset"))
	}
	if (response$status_code != 200) {
		return(webRequestError(response))
	}
	jsonData <- content(response)
	idTypes <- jsonData[[organism]]
	names <- unlist(lapply(idTypes,function(e){return(e$name)}))
	standardIds <- unlist(lapply(idTypes,function(e){return(e$type)}))
	idTypes <- data.frame(name=names, standardId=standardIds, stringsAsFactors=FALSE)
	return(filter(idTypes, name == idType)[[1, 2]])
}


idMappingOutput <- function(mappingOutput,outputFileName,unMapF,dataType,mappingList,sourceIdType,targetIdType){
	if(mappingOutput==TRUE){
		if(length(unMapF)>0){
			write.table(unMapF,file=paste(outputFileName,"_unmappedList.txt",sep=""),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
		}
		if(dataType=="gmt"){
			x <- tapply(mappingList[,targetIdType],mappingList[,1],paste,collapse="\t")
			y <- unique(mappingList[,c(1,2)])
			x <- x[y[,1]]
			z <- data.frame(geneset=y[,1],link=y[,2],gene=x,stringsAsFactors=FALSE)
			write.table(z,file=paste(outputFileName,"_mappedList_from_",sourceIdType,"_to_",targetIdType,".gmt",sep=""),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
		}else{
			if(dataType=="rnk"){
				write.table(mappingList,file=paste(outputFileName,"_mappedList_from_",sourceIdType,"_to_",targetIdType,".rnk",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
			}else{
				write.table(mappingList,file=paste(outputFileName,"_mappedList_from_",sourceIdType,"_to_",targetIdType,".txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
			}
		}
	}
}
