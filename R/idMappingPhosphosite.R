idMappingPhosphosite <- function(organism="hsapiens", dataType="list", inputGeneFile=NULL, inputGene=NULL, sourceIdType, targetIdType, standardId="phosphositeSeq", collapseMethod="mean", mappingOutput=FALSE,  outputFileName="", hostName="http://www.webgestalt.org/") {

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
		x <- tapply(inputGene[,2],inputGene[,1],collapseMethod)
		inputGene <- data.frame(id=names(x),score=as.numeric(x),stringsAsFactors=FALSE)
		inputGeneL <- inputGene[,1]
		colnames(inputGene) <- c(sourceIdType,"score")
	}

	if(dataType=="gmt"){
		colnames(inputGene) <- c("geneset", "link", sourceIdType)
		inputGeneL <- unique(inputGene[,3])
	}

	response <- POST(file.path(hostName, "api", "idmapping"), encode="json",
				body=list(organism=organism, sourcetype=sourceIdType,
				targettype=targetIdType, ids=inputGeneL))

	if (response$status_code != 200) {
		return(webRequestError(response))
	}
	mapRes <- content(response)
	if (mapRes$status == 1) {
		return(webApiError(mapRes))
	}

	mappedIds <- mapRes$mapped
	unmappedIds <- mapRes$unmapped

	if (length(mappedIds) == 0) { return(idMappingError("empty")) }

	names <- c("sourceid", "targetid")
	mappedInputGene <- data.frame(matrix(unlist(lapply(mappedIds, FUN=function(x) { x[names] })), nrow=length(mappedIds), byrow=TRUE), stringsAsFactors=FALSE)
	colnames(mappedInputGene) <- c("userid", targetIdType)

	### Get gene name and symbol in 2nd step, either direct by geneid or mapping to uniprot ambiguously
	if (grepl("Uniprot", sourceIdType, fixed=TRUE) || grepl("Ensembl", sourceIdType, fixed=TRUE) || grepl("Refseq", sourceIdType, fixed=TRUE)) { ##if the sourceIdType is Uniprot, Ensembl or Refseq, directly extract the gene level id####
		mappedInputGene$gene <- unlist(lapply(strsplit(mappedInputGene[, "userid"], "_"), .combineG))
	}else{
		###If the input id type is sequence, we will first map the sequence to uniprot. And then map the uniprot to gene name####
		if (targetIdType == "phosphositeUniprot") {
			mappedInputGene$gene <- unlist(lapply(strsplit(mappedInputGene[, targetIdType], "_"), .combineG))
		} else {
			response <- POST(file.path(hostName, "api", "idmapping"), encode="json",
						body=list(organism=organism, sourcetype=sourceIdType,
						targettype="phosphositeUniprot", ids=inputGeneL))

			if (response$status_code != 200) {
				return(webRequestError(response))
			}
			uniMapRes <- content(response)
			if (uniMapRes$status == 1) {
				return(webApiError(uniMapRes))
			}
			if (length(uniMapRes$mapped) == 0) { return(idMappingError("empty")) }

			names <- c("sourceid", "targetid")
			uniMapRes <- data.frame(matrix(unlist(lapply(uniMapRes$mapped, FUN=function(x) { x[names] })), nrow=length(uniMapRes$mapped), byrow=TRUE), stringsAsFactors=FALSE)
			colnames(uniMapRes) <- c("userid", targetIdType)
			uniMapRes$gene <- unlist(lapply(strsplit(uniMapRes[, targetIdType], "_"), .combineG))
			# Map ID may change nrow due to unmapped ones
			mappedInputGene <- merge(x=mappedInputGene, y=uniMapRes[, c("userid", "gene")], by="userid", all.x=TRUE)
		}
	}

	#####Hard code#######
	if (grepl("Uniprot", sourceIdType, fixed=TRUE) || sourceIdType == "phosphositeSeq") {
		geneType <- "uniprot_swissprot"
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
	mappedInputGene$glink <- paste(outLink,mappedInputGene[, "gene"],sep="")

	########Get gene level information#########
	entrezgeneMapRes <- idMappingGene(organism=organism, dataType="list", inputGene=mappedInputGene[, "gene"], sourceIdType=geneType, targetIdType="entrezgene", mappingOutput=FALSE, hostName=hostName)

	mappedGeneInfo <- entrezgeneMapRes$mapped[, c("userid", "genesymbol", "genename")]
	colnames(mappedGeneInfo) <- c("gene", "genesymbol", "genename")
	#mapFG <- mapFG[,-ncol(mapFG)]

	#mergeA <- merge(x=idTypeGMap, y=mapFG, by.x="gene", by.y="userid", all.x=TRUE)
	#mergeA <- mergeA[,c(idType,"genesymbol","genename","glink")]
	#mergeB <- merge(x=mapF, y=mergeA, by.x="sourceid", by.y=idType, all.x=TRUE)

	mergedRes <- merge(x=mappedInputGene, y=mappedGeneInfo, by="gene", all.x=TRUE)

	if(dataType=="list"){
		inputGene <- mergedRes[, c("userid", "genesymbol", "genename", targetIdType, "glink")]
	}

	if(dataType=="rnk"){
		inputGene <- merge(x=mergedRes, y=inputGene, by.x="userid", by.y=sourceIdType, all.x=TRUE)
		inputGene <- inputGene[,c("userid", "genesymbol", "genename", targetIdType, "score", "glink")]
	}

	if(dataType=="gmt"){
		inputGene <- merge(x=mergedRes, y=inputGene, by.x="userid", by.y=sourceIdType, all.x=TRUE)
		inputGene <- inputGene[, c("geneset", "link", "userid", "genesymbol", "genename", targetIdType, "glink")]
	}

	#############Output#######################
	idMappingOutput(mappingOutput, outputFileName, unmappedIds, dataType, inputGene, sourceIdType, targetIdType=targetIdType)
	r <- list(mapped=inputGene,unmapped=unmappedIds)
	return(r)
}


.combineG <- function(e){
	e <- e[-length(e)]
	e <- paste(e,collapse="_")
	return(e)
}
