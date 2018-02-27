IDMapping_gene <- function(organism="hsapiens",dataType="list",inputGeneFile=NULL,inputGene=NULL,sourceIdType,standardID,targetIdType,collapseMethod="mean",mappingOutput=FALSE, outputFileName="",methodType="R",hostName="http://www.webgestalt.org/"){
	goldIdType <- c("entrezgene","genesymbol","genename")

	largeIdList <- fread(input=file.path(hostName,"data","largeIdList.txt"),header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses="character",data.table=FALSE,showProgress=FALSE)
	largeIdList <- as.character(largeIdList[,1])

	###########Check input data type###############

	inputGene <- IDMapping_input(dataType=dataType,inputGeneFile=inputGeneFile,inputGene=inputGene)
	if(.hasError(inputGene)){
		return(inputGene)
	}

	##########ID Mapping Specify to gene level###############

	re <- .processSourceIDMapGene(hostName=hostName,organism=organism,largeIdList=largeIdList,inputGene=inputGene,standardID=standardID,dataType=dataType,idType=sourceIdType,collapseMethod=collapseMethod,methodType=methodType)
	if(.hasError(re)){
		return(re)
	}
	inputGene <- re$mapped
	unMapF <- re$unmapped

	###########If the target ID type is not from goldIdType#######
	inputGene <- unique(inputGene)

	if(length(which(goldIdType==targetIdType))==0 && sourceIdType!=targetIdType){
		x <- unique(inputGene[,standardID])
		targetF <- IDMapping_map(largeIdList=largeIdList,sourceIdType=targetIdType,standardID=standardID,hostName=hostName,organism=organism,inputGene=x,mapType="target")
		if(.hasError(targetF)){
			return(targetF)
		}
		targetF <- targetF$mapF

		inputGene1 <- merge(x=inputGene,y=targetF,by=standardID)
		inputGene1 <- unique(inputGene1)
		unMapF2 <- setdiff(inputGene[,"userid"],inputGene1[,"userid"])
		unMapF <- unique(c(unMapF,unMapF2))
		inputGene <- inputGene1
		if(dataType=="list"){
			inputGene <- inputGene[,c(2,3,4,1,6,5)]
			colnames(inputGene) <- c("userid","genesymbol","genename","entrezgene",targetIdType,"glink")
		}
		if(dataType=="rnk"){
			inputGene <- inputGene[,c(2,3,4,1,7,5,6)]
			colnames(inputGene) <- c("userid","genesymbol","genename","entrezgene",targetIdType,"score","glink")
		}
		if(dataType=="gmt"){
			inputGene <- inputGene[,c(2:6,1,8,7)]
			colnames(inputGene) <- c("geneset","link","userid","genesymbol","genename","entrezgene",targetIdType,"glink")
		}
	}

	#############Output#######################

	IDMapping_output(mappingOutput,outputFileName,unMapF,dataType,inputGene,sourceIdType,targetIdType)
	r <- list(mapped=inputGene,unmapped=unMapF)
	return(r)
}


.processSourceIDMapGene <- function(hostName,organism,largeIdList,inputGene,standardID,dataType,idType,collapseMethod,methodType){
	geneSymbol <- fread(input=file.path(hostName,"data","xref",paste(organism,"_genesymbol_",standardID,".table",sep="")),header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses="character",data.table=FALSE,showProgress=FALSE)

	geneName <- fread(input=file.path(hostName,"data","xref",paste(organism,"_genename_",standardID,".table",sep="")),header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses="character",data.table=FALSE,showProgress=FALSE)

	colnames(geneSymbol) <- c("entrezgeneS","genesymbol")

	colnames(geneName) <- c("entrezgeneS","genename")

	geneAnn <- merge(x=geneSymbol,y=geneName,by="entrezgeneS",all=TRUE)

	if(dataType=="list"){
		inputGeneL <- unique(inputGene)
	}

	if(dataType=="rnk"){
		######Collapse the gene ids with multiple scores##########
		x <- tapply(inputGene[,2],inputGene[,1],collapseMethod)
		inputGene <- data.frame(id=names(x),score=as.numeric(x),stringsAsFactors=FALSE)
		inputGeneL <- inputGene[,1]
		colnames(inputGene) <- c(idType,"score")
	}

	if(dataType=="gmt"){
		inputGeneL <- unique(inputGene[,3])
	}

	mapR <- IDMapping_map(largeIdList=largeIdList,sourceIdType=idType,standardID=standardID,hostName=hostName,organism=organism,inputGene=inputGeneL,mapType="source",methodType=methodType)
	if(.hasError(mapR)){
		return(mapR)
	}

	mapF <- mapR$mapF
	unMapF <- mapR$unMapF

	colnames(mapF) <- c("entrezgeneS",idType)

	if(dataType=="list"){
		inputGene <- merge(x=mapF,y=geneAnn,by="entrezgeneS",all.x=TRUE)
		inputGene <- inputGene[,c(2,3,4,1)]
		colnames(inputGene) <- c("userid","genesymbol","genename","entrezgene")
	}

	if(dataType=="rnk"){
		colnames(inputGene) <- c(idType,"score")
		inputGene <- merge(x=mapF,y=inputGene,by=idType,all.x=TRUE)
		inputGene <- merge(x=inputGene,y=geneAnn,by="entrezgeneS",all.x=TRUE)
		inputGene <- inputGene[,c(2,4,5,1,3)]
		colnames(inputGene) <- c("userid","genesymbol","genename","entrezgene","score")
	}

	if(dataType=="gmt"){
		colnames(inputGene) <- c("geneset","link",idType)
		inputGene <- merge(x=mapF,y=inputGene,by=idType,all.x=TRUE)
		inputGene <- merge(x=inputGene,y=geneAnn,by="entrezgeneS",all.x=TRUE)
		inputGene <- as.matrix(inputGene)
		inputGene <- inputGene[,c(3,4,2,5,6,1)]
		colnames(inputGene) <- c("geneset","link","userid","genesymbol","genename","entrezgene")
	}

	inputGene <- as.data.frame(inputGene)
	inputGene$glink <- paste("https://www.ncbi.nlm.nih.gov/gene/?term=",inputGene[,"entrezgene"],sep="")

	re <- list(mapped=inputGene,unmapped=unMapF)
	return(re)
}
