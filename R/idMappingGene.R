idMappingGene <- function(organism="hsapiens", dataType="list", inputGeneFile=NULL, inputGene=NULL, sourceIdType, targetIdType, collapseMethod="mean", mappingOutput=FALSE,  outputFileName="", hostName="http://www.webgestalt.org/") {

	###########Check input data type###############
	inputGene <- idMappingInput(dataType=dataType,inputGeneFile=inputGeneFile,inputGene=inputGene)
	if(.hasError(inputGene)){
		return(inputGene)
	}

	##########ID Mapping Specify to gene level###############
	if(dataType=="list"){
		inputGeneL <- unique(inputGene)
	}

	if(dataType=="rnk"){
		######Collapse the gene ids with multiple scores##########
		x <- tapply(inputGene[,2],inputGene[,1],collapseMethod)
		inputGene <- data.frame(id=names(x),score=as.numeric(x),stringsAsFactors=FALSE)
		inputGeneL <- inputGene[,1]
		colnames(inputGene) <- c(sourceIdType,"score")
	}

	if(dataType=="gmt"){
		colnames(inputGene) <- c("geneset", "link", sourceIdType)
		inputGeneL <- unique(inputGene[,3])
	}

	mapR <- POST(file.path(hostName, "api", "idmapping"), encode="json",
				body=list(organism=organism, sourcetype=sourceIdType,
				targettype=targetIdType, ids=inputGeneL))

	if (mapR$status_code != 200) {
		return(webRequestError(mapR))
	}
	mapR <- content(mapR)
	if (mapR$status == 1) {
		return(webApiError(mapR))
	}

	if(.hasError(mapR)){
		return(mapR)
	}

	mappedIds <- mapR$mapped
	unmappedIds <- mapR$unmapped

	if (length(mappedIds) == 0) { return(idMappingError("empty")) }

	names <- c("sourceid", "genesymbol", "genename", "targetid")
	mappedInputGene <- data.frame(matrix(unlist(lapply(mappedIds, FUN=function(x) { x[names] })), nrow=length(mappedIds), byrow=TRUE), stringsAsFactors=FALSE)
	colnames(mappedInputGene) <- c("userid", "genesymbol","genename", targetIdType)
	if (dataType=="list") {
		inputGene <- mappedInputGene
	} else if (dataType=="rnk") {
		inputGene <- merge(x=mappedInputGene, y=inputGene, by.x="userid", by.y=sourceIdType)
	} else if (dataType=="gmt") {
		inputGene <- merge(x=mappedInputGene, y=inputGene, by.x="userid", by.y=sourceIdType)
		inputGene <- inputGene[, c("geneset", "link", "userid", "genesymbol", "genename", targetIdType)]
	}

	if (targetIdType != "entrezgene" && sourceIdType!=targetIdType) {
		entrezgeneMapRes <- idMappingGene(organism, dataType="list", inputGene=inputGeneL, sourceIdType=sourceIdType, targetIdType="entrezgene")
		entrezgeneMapInfp <- entrezgeneMapRes$mapped[, c("userid", "entrezgene")]
		inputGene <- merge(x=inputGene, y=entrezgeneMapRes, by="userid", all.x=TRUE)
		if (dataType=="list") {
			inputGene <- inputGene[, c("userid", "genesymbol", "genename", "entrezgene", targetIdType)]
		} else if (dataType=="rnk") {
			inputGene <- inputGene[, c("userid", "genesymbol", "genename", "entrezgene", targetIdType, "score")]
		} else if (dataType=="gmt") {
			inputGene <- inputGene[, c("geneset", "link", "userid", "genesymbol", "genename", "entrezgene", targetIdType)]
		}
	}
	inputGene$glink <- paste0("https://www.ncbi.nlm.nih.gov/gene/?term=", inputGene[, "entrezgene"])

	#############Output#######################
	idMappingOutput(mappingOutput,outputFileName,unMapF,dataType,inputGene,sourceIdType,targetIdType)
	r <- list(mapped=inputGene,unmapped=unmappedIds)
	return(r)
}
