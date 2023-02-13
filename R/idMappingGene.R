#' @importFrom httr POST content
#' @importFrom dplyr inner_join select filter left_join %>%
idMappingGene <- function(organism="hsapiens", dataType="list", inputGeneFile=NULL, inputGene=NULL, sourceIdType, targetIdType, collapseMethod="mean", mappingOutput=FALSE,  outputFileName="", hostName="https://www.webgestalt.org/") {

	###########Check input data type###############
	inputGene <- idMappingInput(dataType=dataType,inputGeneFile=inputGeneFile,inputGene=inputGene)

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
	if (startsWith(hostName, "file://")) {
		# old way of mapping with mapping files. Now only used for WebGestaltReporter when hostName is file protocol
		sourceMap <- read_tsv(
			removeFileProtocol(file.path(hostName, "xref", paste(organism, sourceIdType, "entrezgene.table", sep="_"))),
			col_names=c("entrezgene", "userId"), col_types="cc", quote=""
		) %>% filter(.data$userId %in% inputGeneL)
		symbolMap <- read_tsv(
			removeFileProtocol(file.path(hostName, "xref", paste(organism, "genesymbol", "entrezgene.table", sep="_"))),
			col_names=c("entrezgene", "geneSymbol"), col_types="cc", quote=""
		)
		nameMap <- read_tsv(
			removeFileProtocol(file.path(hostName, "xref", paste(organism, "genename", "entrezgene.table", sep="_"))),
			col_names=c("entrezgene", "geneName"), col_types="cc", quote=""
		)
		sourceMap <- sourceMap %>% left_join(symbolMap, by=c("entrezgene")) %>% left_join(nameMap, by=c("entrezgene"))

		if (targetIdType %in% c("entrezgene", sourceIdType)) {
			mappedInputGene <- sourceMap
		} else {
			targetMap <- read_tsv(removeFileProtocol(file.path(hostName, "xref", paste(organism, targetIdType, "entrezgene.table", sep="_"))),
				col_names=c("entrezgene", targetIdType), col_types="cc", quote="")
			mappedInputGene <- inner_join(sourceMap, targetMap, by=c("entrezgene"))
		}
		if (nrow(mappedInputGene) == 0) { return(idMappingError("empty")) }
		mappedInputGene <- mappedInputGene %>%
			select(.data$userId, .data$geneSymbol, .data$geneName, targetIdType)
		unmappedIds <- setdiff(inputGeneL, mappedInputGene$userId)

	} else {
		# new way uses web server API
		mapR <- POST(file.path(hostName, "api", "idmapping"), encode="json",
			body=list(organism=organism, sourceType=sourceIdType,
			targetType=targetIdType, ids=inputGeneL)
		)

		if (mapR$status_code != 200) {
			stop(webRequestError(mapR))
		}
		mapR <- content(mapR)
		if (mapR$status == 1) {
			stop(webApiError(mapR))
		}

		mappedIds <- mapR$mapped
		unmappedIds <- unlist(mapR$unmapped)
		if (is.null(targetIdType)) {
			targetIdType <- mapR$standardId
		}

		if (length(mappedIds) == 0) { stop(idMappingError("empty")) }

		names <- c("sourceId", "geneSymbol", "geneName", "targetId")
		mappedInputGene <- data.frame(matrix(unlist(lapply(replace_null(mappedIds), FUN=function(x) { x[names] })), nrow=length(mappedIds), byrow=TRUE), stringsAsFactors=FALSE)
		colnames(mappedInputGene) <- c("userId", "geneSymbol", "geneName", targetIdType)
	}

	if (dataType=="list") {
		inputGene <- mappedInputGene
	} else if (dataType=="rnk") {
		inputGene <- inner_join(mappedInputGene, inputGene, by=c("userId"=sourceIdType))
	} else if (dataType=="gmt") {
		inputGene <- inner_join(mappedInputGene, inputGene, by=c("userId"=sourceIdType)) %>%
			select(.data$geneSet, .data$description, .data$userId, .data$geneSymbol, .data$geneName, targetIdType)
	}

	if (targetIdType != "entrezgene" && sourceIdType!=targetIdType) {
		entrezgeneMapRes <- idMappingGene(organism, dataType="list", inputGene=inputGeneL, sourceIdType=sourceIdType, targetIdType="entrezgene", hostName=hostName)
		inputGene <- left_join(inputGene, entrezgeneMapRes$mapped, by="userId")

		if (dataType=="list") {
			inputGene <- select(inputGene, .data$userId, geneSymbol=.data$geneSymbol.x, geneName=.data$geneName.x, .data$entrezgene, targetIdType)
		} else if (dataType=="rnk") {
			inputGene <- select(inputGene, .data$userId, geneSymbol=.data$geneSymbol.x, geneName=.data$geneName.x, .data$entrezgene, targetIdType, .data$score)
		} else if (dataType=="gmt") {
			inputGene <- select(inputGene, .data$geneSet, .data$description, .data$userId, geneSymbol=.data$geneSymbol.x, geneName=.data$geneName.x, .data$entrezgene, targetIdType)
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
