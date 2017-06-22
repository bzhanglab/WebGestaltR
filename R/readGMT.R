readGMT <- function(gmtFile){
#####Change a gmt file to a three column matrix (gene set name, gene set description and genes)#######
	if(file_extension(gmtFile)!="gmt"){
		return(gmtFormatError("empty"))
	}else{
		data <- readLines(gmtFile)
		data <- strsplit(data,"\t")
		data <- lapply(data,.toList)
		data <- do.call("rbind",data)
		
		if(is.null(data)){
			return(gmtFormatError("incorrect"))
		}else{
			return(data)
		}
	}
}


.toList <- function(data){
	if(length(data)>2){
		data <- data[!is.na(data)]
		data1 <- cbind(rep(data[1],length(data)-2),rep(data[2],length(data)-2),data[c(-1,-2)])
		return(data1)
	}else{
		return(NULL)
	}
}