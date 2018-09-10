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
		colnames(inputGene) <- c("geneset", "link", sourceIdType)
		inputGeneL <- unique(inputGene$gene)
	}

	response <- POST(file.path(hostName, "api", "idmapping"), encode="json",
				body=list(organism=organism, sourcetype=sourceIdType,
				targettype=targetIdType, ids=inputGeneL, standardid="phosphositeSeq"
				))

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

	names <- c("sourceid", "targetid")
	mappedInputGene <- data.frame(matrix(unlist(lapply(mappedIds, FUN=function(x) { x[names] })), nrow=length(mappedIds), byrow=TRUE), stringsAsFactors=FALSE)
	colnames(mappedInputGene) <- c("userid", targetIdType)

	### Get gene name and symbol in 2nd step, either direct by geneid or mapping to uniprot ambiguously
	if (grepl("Uniprot", sourceIdType, fixed=TRUE) || grepl("Ensembl", sourceIdType, fixed=TRUE) || grepl("Refseq", sourceIdType, fixed=TRUE)) { ##if the sourceIdType is Uniprot, Ensembl or Refseq, directly extract the gene level id####
		mappedInputGene$gene <- unlist(lapply(strsplit(mappedInputGene$userid, "_"), .combineG))
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
			mappedInputGene <- uniMapRes %>% select(userid, gene) %>% right_join(mappedInputGene, by=userid)
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
	mappedInputGene$glink <- paste0(outLink, mappedInputGene$gene)

	########Get gene level information#########
	entrezgeneMapRes <- idMappingGene(organism=organism, dataType="list", inputGene=mappedInputGene$gene, sourceIdType=geneType, targetIdType="entrezgene", mappingOutput=FALSE, hostName=hostName)

	mergedRes <- entrezgeneMapRes$mapped %>% select(gene=userid, genesymbol, genename) %>%
		right_join(mappedInputGene, by="gene")

	if(dataType=="list"){
		inputGene <- select(mergedRes, userid, genesymbol, genename, targetIdType, glink)
	}

	if(dataType=="rnk"){
		inputGene <- mergedRes %>% left_join(inputGene, by=c("userid"=sourceIdType)) %>%
			select(userid, genesymbol, genename, targetIdType, score, glink)
	}

	if(dataType=="gmt"){
		inputGene <- mergedRes %>% left_join(inputGene, by=c("userid"=sourceIdType)) %>%
			select(geneset, link, userid, genesymbol, genename, targetIdType, glink)
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
