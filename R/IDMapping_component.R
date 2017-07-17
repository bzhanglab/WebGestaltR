IDMapping_input <- function(dataType="list",inputGeneFile,inputGene){
	 if(dataType=="gmt"){
    	if(!is.null(inputGeneFile)){
    		inputGene <- readGMT(inputGeneFile)
        return(inputGene)
    	}else{
    		return(gmtFormatError("empty"))
    	}
    }else{
    	inputGene <- formatCheck(dataType=dataType,inputGeneFile=inputGeneFile,inputGene=inputGene)
    	return(inputGene)
    }	
}


identifyStandardId <- function(hostName,idtype,organism,type){
		if(type=="interest"){
			json_data <- fromJSON(file=file.path(hostName,"data","idtypesummary.json"))
		}
		if(type=="reference"){
			json_data <- fromJSON(file=file.path(hostName,"data","referenceSetsummary.json"))
		}
		idtypes <- json_data[[organism]]
		names <- unlist(lapply(idtypes,function(e){return(e$name)}))
		staIDs <- unlist(lapply(idtypes,function(e){return(e$type)}))
		idtypes <- cbind(names,staIDs)
		return(idtypes[idtypes[,1]==idtype,2])
}


IDMapping_map <- function(largeIdList,sourceIdType,standardID,hostName,organism,inputGene,mapType,methodType){
	if(length(which(largeIdList==sourceIdType))>0){
		if(methodType=="Python"){
	  	re <- .Python_IDMap(hostName,organism,sourceIdType,standardID,inputGene,mapType)
      if(.hasError(re)){
	  		return(re)
	  	}
	  }else{
	  	re <- .R_IDMap(hostName,organism,sourceIdType,standardID,inputGene,mapType)
	  }
	}else{
		re <- .R_IDMap(hostName,organism,sourceIdType,standardID,inputGene,mapType)
	}
	return(re)
}

.R_IDMap <- function(hostName,organism,sourceIdType,standardID,inputGene,mapType){
#if mapType is source, inputGene is other ids. If mapType is target, inputGene is standardID

		  sourceFile <- fread(input=file.path(hostName,"data","xref",paste(organism,"_",sourceIdType,"_",standardID,".table",sep="")),header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses="character",data.table=FALSE)
		  colnames(sourceFile) <- c(standardID,sourceIdType)
		  topF <- paste(sourceFile[1:5,2],collapse=",")
		  if(mapType=="source"){
				mapF <- sourceFile[sourceFile[,2] %in% inputGene,,drop=FALSE]
				unMapF <- setdiff(inputGene,mapF[,2])
			}else{
				mapF <- sourceFile[sourceFile[,1] %in% inputGene,,drop=FALSE]
				unMapF <- setdiff(inputGene,mapF[,1])
			}
			
			if(nrow(mapF)==0){
					return(IDMappingError(type="unmapped",idType=sourceIdType,topF=topF))
			}else{
				re <- list(mapF=mapF,unMapF=unMapF)
				return(re)
			}
}

.Python_IDMap <- function(hostName,organism,sourceIdType,standardID,inputGene,mapType){

#for mac, install pip first sudo easy_install pip
#install pandas module, 
	if(pyIsConnected()){
			re <- tryCatch(pyExec('import pandas as pd'),error=function(e){return(NULL)})
			if(is.null(re)){
				return(IDMappingError(type="pandas"))
			}
			pySet('hostName',hostName)
			pySet('organism',organism)
			pySet('sourceIdType',sourceIdType)
			pySet('standardID',standardID)
			pySet('inputGene',as.list(inputGene))
			
			pyExec('data = pd.read_csv(filepath_or_buffer=hostName+"/data/xref/"+organism+"_"+sourceIdType+"_"+standardID+".table",sep="\t",memory_map=True,header=None,names=[standardID,sourceIdType],dtype=str)')
			pyExec('top = data.head(n=5)')
			pyExec('top = top.to_dict()')
			topF <- PythonInR::pyGet('top')
			topF <-do.call(cbind,topF)
			topF <- topF[,c(standardID,sourceIdType)]
			topF <- topF[,2]
			topF <- paste(topF,collapse=",")
			if(mapType=="source"){
				pyExec('map = data.loc[data[sourceIdType].isin(inputGene)]')
			}else{
				pyExec('map = data.loc[data[standardID].isin(inputGene)]')
			}
			pyExec('map = map.to_dict()')
			mapF <- pyGet('map')
			mapF <- do.call(cbind,mapF)
			mapF <- mapF[,c(standardID,sourceIdType)]
			
			if(is.null(mapF)){
					return(IDMappingError(type="unmapped",idType=sourceIdType,topF=topF))
			}else{
				if(mapType=="source"){
					unMapF <- setdiff(inputGene,mapF[,2])
				}else{
					unMapF <- setdiff(inputGene,mapF[,1])
				}
				re <- list(mapF=mapF,unMapF=unMapF)
				return(re)
			}
	}else{
		return(IDMappingError(type="python"))
	}
}

IDMapping_output <- function(mappingOutput,outputFileName,unMapF,dataType,mappingList,sourceIdType,targetIdType){
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
