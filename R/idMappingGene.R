#' @importFrom httr POST content
#' @importFrom dplyr inner_join select left_join %>%
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
		x <- tapply(inputGene$score, inputGene$gene, collapseMethod)
		inputGene <- data.frame(gene=names(x),score=as.numeric(x),stringsAsFactors=FALSE)
		inputGeneL <- inputGene$gene
		colnames(inputGene) <- c(sourceIdType,"score")
	}

	if(dataType=="gmt"){
		colnames(inputGene) <- c("geneSet", "description", sourceIdType)
		inputGeneL <- unique(inputGene[[sourceIdType]])
	}

	mapR <- POST(file.path(hostName, "api", "idmapping"), encode="json",
				body=list(organism=organism, sourceType=sourceIdType,
				targetType=targetIdType, ids=inputGeneL))

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
	unmappedIds <- unlist(mapR$unmapped)
	if (is.null(targetIdType)) {
		targetIdType <- mapR$standardId
	}

	if (length(mappedIds) == 0) { return(idMappingError("empty")) }

	names <- c("sourceId", "geneSymbol", "geneName", "targetId")
	mappedInputGene <- data.frame(matrix(unlist(lapply(mappedIds, FUN=function(x) { x[names] })), nrow=length(mappedIds), byrow=TRUE), stringsAsFactors=FALSE)
	colnames(mappedInputGene) <- c("userId", "geneSymbol", "geneName", targetIdType)
	if (dataType=="list") {
		inputGene <- mappedInputGene
	} else if (dataType=="rnk") {
		inputGene <- inner_join(mappedInputGene, inputGene, by=c("userId"=sourceIdType))
	} else if (dataType=="gmt") {
		inputGene <- inner_join(mappedInputGene, inputGene, by=c("userId"=sourceIdType)) %>%
			select(.data$geneSet, .data$description, .data$userId, .data$geneSymbol, .data$geneName, .data$targetIdType)
	}

	if (targetIdType != "entrezgene" && sourceIdType!=targetIdType) {
		entrezgeneMapRes <- idMappingGene(organism, dataType="list", inputGene=inputGeneL, sourceIdType=sourceIdType, targetIdType="entrezgene", hostName=hostName)
		inputGene <- left_join(inputGene, entrezgeneMapRes$mapped, by="userId")

		if (dataType=="list") {
			inputGene <- select(inputGene, .data$userId, geneSymbol=.data$geneSymbol.x, geneName=.data$geneName.x, .data$entrezgene, .data$targetIdType)
		} else if (dataType=="rnk") {
			inputGene <- select(inputGene, .data$userId, geneSymbol=.data$geneSymbol.x, geneName=.data$geneName.x, .data$entrezgene, .data$targetIdType, .data$score)
		} else if (dataType=="gmt") {
			inputGene <- select(inputGene, .data$geneSet, .data$description, .data$userId, geneSymbol=.data$geneSymbol.x, geneName=.data$geneName.x, .data$entrezgene, .data$targetIdType)
		}
	}
	inputGene$gLink <- paste0("https://www.ncbi.nlm.nih.gov/gene/?term=", inputGene$entrezgene)

	#############Output#######################
	if (mappingOutput) {
		idMappingOutput(outputFileName, inputGene, unmappedIds, dataType, sourceIdType, targetIdType)
	}
	r <- list(mapped=inputGene,unmapped=unmappedIds)
	return(r)
}
