#' @importFrom httr POST content
#' @importFrom dplyr right_join select left_join %>%
idMappingPhosphosite <- function(organism="hsapiens", dataType="list", inputGeneFile=NULL, inputGene=NULL, sourceIdType, targetIdType, collapseMethod="mean", mappingOutput=FALSE,  outputFileName="", hostName="http://www.webgestalt.org/") {

	###########Check input data type###############
	inputGene <- idMappingInput(dataType=dataType,inputGeneFile=inputGeneFile,inputGene=inputGene)
	if(.hasError(inputGene)){
		return(inputGene)
	}

	##########ID Mapping Specify to phosphosite level###############
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
		colnames(inputGene) <- c("geneSet", "link", sourceIdType)
		inputGeneL <- unique(inputGene$gene)
	}

	if (startsWith(hostName, "file://")) {
		sourceMap <- read_tsv(
			removeFileProtocol(file.path(hostName, "xref", paste(organism, sourceIdType, "phosphositeSeq.table", sep="_"))),
			col_names=c("phosphositeSeq", "userId"), col_types="cc", quote=""
		) %>% filter(.data$userId %in% inputGeneL)
		if (targetIdType == "phosphositeSeq" || targetIdType == sourceIdType) {
			mappedInputGene <- sourceMap
		} else {
			targetMap <- read_tsv(
				removeFileProtocol(file.path(hostName, "xref", paste(organism, targetIdType, "phosphositeSeq.table", sep="_"))),
				col_names=c("phosphositeSeq", targetIdType), col_types="cc", quote=""
			)
			mappedInputGene <- inner_join(sourceMap, targetMap, by=c("phosphositeSeq"))
		}
		if (nrow(mappedInputGene) == 0) { return(idMappingError("empty")) }
		mappedInputGene <- select(mappedInputGene, .data$userId, targetIdType)
		unmappedIds <- setdiff(inputGeneL, mappedInputGene$userId)
	} else {
		response <- POST(file.path(hostName, "api", "idmapping"), encode="json",
			body=list(organism=organism, sourceType=sourceIdType,
			targetType=targetIdType, ids=inputGeneL, standardId="phosphositeSeq")
		)
		if (response$status_code != 200) {
			return(webRequestError(response))
		}
		mapRes <- content(response)
		if (mapRes$status == 1) {
			return(webApiError(mapRes))
		}

		mappedIds <- mapRes$mapped
		unmappedIds <- unlist(mapRes$unmapped)

		if (length(mappedIds) == 0) { return(idMappingError("empty")) }
		names <- c("sourceId", "targetId")
		mappedInputGene <- data.frame(matrix(unlist(lapply(mappedIds, FUN=function(x) { x[names] })), nrow=length(mappedIds), byrow=TRUE), stringsAsFactors=FALSE)
		colnames(mappedInputGene) <- c("userId", targetIdType)
	}

	### Get gene name and symbol in 2nd step, either direct by geneid or mapping to uniprot ambiguously
	# TODO mapping to target other than 15mer may introduce ambiguity, like DTQIKRNtFVGTPFW maps to three STKs in uniprot.
	# Not essential for WG, but could use protein ID to determine
	if (grepl("Uniprot", sourceIdType, fixed=TRUE) || grepl("Ensembl", sourceIdType, fixed=TRUE) || grepl("Refseq", sourceIdType, fixed=TRUE)) { ##if the sourceIdType is Uniprot, Ensembl or Refseq, directly extract the gene level id####
		mappedInputGene$gene <- unlist(lapply(strsplit(mappedInputGene$userId, "_"), .combineG))
	}else{
		###If the input id type is sequence, we will first map the sequence to uniprot. And then map the uniprot to gene name####
		if (targetIdType == "phosphositeUniprot") {
			mappedInputGene$gene <- unlist(lapply(strsplit(mappedInputGene[, targetIdType], "_"), .combineG))
		} else {
			if (startsWith(hostName, "file://")) {
				uniMapRes <- read_tsv(
					removeFileProtocol(file.path(hostName, "xref", paste(organism, "phosphositeUniprot", "phosphositeSeq.table", sep="_"))),
					col_names=c("phosphositeSeq", "phosphositeUniprot"), col_types="cc", quote=""
				) %>% filter(.data$phosphositeSeq %in% mappedInputGene$phosphositeSeq)
			} else {
				response <- POST(file.path(hostName, "api", "idmapping"), encode="json",
					body=list(organism=organism, sourceType="phosphositeSeq", standardId="phosphositeSeq",
					targetType="phosphositeUniprot", ids=inputGeneL)
				)

				if (response$status_code != 200) {
					return(webRequestError(response))
				}
				uniMapRes <- content(response)
				if (uniMapRes$status == 1) {
					return(webApiError(uniMapRes))
				}
				if (length(uniMapRes$mapped) == 0) { return(idMappingError("empty")) }

				names <- c("sourceId", "targetId")
				uniMapRes <- data.frame(matrix(unlist(lapply(uniMapRes$mapped, FUN=function(x) { x[names] })), nrow=length(uniMapRes$mapped), byrow=TRUE), stringsAsFactors=FALSE)
				colnames(uniMapRes) <- c("phosphositeSeq", "phosphositeUniprot")
			}

			uniMapRes$gene <- unlist(lapply(strsplit(uniMapRes[, "phosphositeUniprot"], "_"), .combineG))
			# Map ID may change nrow due to unmapped ones
			mappedInputGene <- uniMapRes %>% select(.data$phosphositeSeq, .data$gene) %>% right_join(mappedInputGene, by="phosphositeSeq")
		}
	}
	#####Hard code#######
	if (grepl("Uniprot", sourceIdType, fixed=TRUE) || sourceIdType == "phosphositeSeq") {
		geneType <- "uniprotswissprot"
		outLink <- "http://www.uniprot.org/uniprot/"
	}

	if (grepl("Ensembl", sourceIdType, fixed=TRUE)) {
		geneType <- "ensembl_peptide_id"
		outLink <- paste("http://www.ensembl.org/",organism,"/Gene/Summary?db=core;t=",sep="")
	}

	if (grepl("Refseq", sourceIdType, fixed=TRUE)) {
		geneType <- "refseq_peptide"
		outLink <- "https://www.ncbi.nlm.nih.gov/protein/"
	}
	mappedInputGene$gLink <- paste0(outLink, mappedInputGene$gene)

	########Get gene level information#########
	entrezgeneMapRes <- idMappingGene(organism=organism, dataType="list", inputGene=mappedInputGene$gene, sourceIdType=geneType, targetIdType="entrezgene", mappingOutput=FALSE, hostName=hostName)

	mergedRes <- entrezgeneMapRes$mapped %>% select(gene=.data$userId, .data$geneSymbol, .data$geneName) %>%
		right_join(mappedInputGene, by="gene")

	if(dataType=="list"){
		inputGene <- select(mergedRes, .data$userId, .data$geneSymbol, .data$geneName, targetIdType, .data$gLink)
	}

	if(dataType=="rnk"){
		inputGene <- mergedRes %>% left_join(inputGene, by=c("userId"=sourceIdType)) %>%
			select(.data$userId, .data$geneSymbol, .data$geneName, targetIdType, .data$score, .data$gLink)
	}

	if(dataType=="gmt"){
		inputGene <- mergedRes %>% left_join(inputGene, by=c("userId"=sourceIdType)) %>%
			select(.data$geneSet, .data$link, .data$userId, .data$geneSymbol, .data$geneName, targetIdType, .data$gLink)
	}

	#############Output#######################
	if (mappingOutput) {
		idMappingOutput(outputFileName, inputGene, unmappedIds, dataType, sourceIdType, targetIdType=targetIdType)
	}
	r <- list(mapped=inputGene,unmapped=unmappedIds)
	return(r)
}


.combineG <- function(e){
	e <- e[-length(e)]
	e <- paste(e,collapse="_")
	return(e)
}
