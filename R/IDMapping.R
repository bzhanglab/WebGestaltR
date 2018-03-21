idMapping <- function(organism="hsapiens",dataType="list",inputGeneFile=NULL,inputGene=NULL,sourceIdType,targetIdType=NULL,collapseMethod="mean",mappingOutput=FALSE, outputFileName="",methodType="R",hostName="http://www.webgestalt.org/"){
	#############Check general parameters########
	errorTest <- parameterErrorMessage(organism=organism,dataType=dataType,collapseMethod=collapseMethod,methodType=methodType,hostName=hostName,mappingOutput=mappingOutput)

	if(!is.null(errorTest)){
		return(errorTest)
	}

	############Check source id type#########
	errorTest <- idTypeError(idType=sourceIdType,organism=organism,hostName=hostName)
	if(!is.null(errorTest)){
		return(errorTest)
	}

	##########Identify the standardId for the input ID type###########
	standardSource <- identifyStandardId(hostName=hostName,idType=sourceIdType,organism=organism,type="interest")

	############Check target id type#########
	if(!is.null(targetIdType)){
		errorTest <- targetIdTypeError(idType=targetIdType,organism=organism,hostName=hostName)
		if(!is.null(errorTest)){
			return(errorTest)
		}else{
			standardTarget <- identifyStandardId(hostName=hostName,idType=targetIdType,organism=organism,type="interest")
			errorTest <- stardardDiffError(standardSource=standardSource,standardTarget=standardTarget)
			if(!is.null(errorTest)){
				return(errorTest)
			}
		}
	}else{
		targetIdType <- standardSource
	}

	##########gene level ID Mapping##########
	if(standardSource=="entrezgene"){
		idMap <- idMappingGene(organism=organism,dataType=dataType,inputGeneFile=inputGeneFile,inputGene=inputGene,sourceIdType=sourceIdType,standardId=standardSource,targetIdType=targetIdType,collapseMethod=collapseMethod,mappingOutput=mappingOutput,outputFileName=outputFileName,methodType=methodType,hostName=hostName)
	}

	if(standardSource=="phosphositeSeq"){
		idMap <- idMappingPhosphosite(organism=organism,dataType=dataType,inputGeneFile=inputGeneFile,inputGene=inputGene,sourceIdType=sourceIdType,standardId=standardSource,collapseMethod=collapseMethod,mappingOutput=mappingOutput,outputFileName=outputFileName,methodType=methodType,hostName=hostName)
	}

	if(.hasError(idMap)){
		return(idMap)
	}else{
		idMap$standardId <- standardSource
		return(idMap)
	}
}

IDMapping <- function(...) {
	cat("WARNING: Function IDMapping is deprecated and changed to idMapping!\n")
	return(idMapping(...))
}
