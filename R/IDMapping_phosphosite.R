IDMapping_phosphosite <- function(organism="hsapiens",dataType="list",inputGeneFile=NULL,inputGene=NULL,sourceIdType,standardID,collapseMethod="mean",mappingOutput=FALSE, outputFileName="",methodType="R",hostName="http://www.webgestalt.org/"){
     
    largeIdList <- fread(input=file.path(hostName,"data","largeIdList.txt"),header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses="character",data.table=FALSE,showProgress=FALSE)
    largeIdList <- as.character(largeIdList[,1])
    
    ###########Check input data type###############
    
    inputGene <- IDMapping_input(dataType=dataType,inputGeneFile=inputGeneFile,inputGene=inputGene)
    if(is.character(inputGene) && length(inputGene)==1 && length(grep("ERROR:",inputGene))>0){
       return(inputGene)
    }
    
    
	 ##########ID Mapping Specify to phosphosite level###############
	 
	  re <- .processSourceIDMapPhosphosite(hostName=hostName,organism=organism,largeIdList=largeIdList,inputGene=inputGene,standardID=standardID,dataType=dataType,idType=sourceIdType,collapseMethod=collapseMethod,methodType=methodType)
		if(is.character(re) && length(re)==1 && length(grep("ERROR:",re))>0){
			return(re)
		}
		inputGene <- re$mapped
		unMapF <- re$unmapped
	  
	  #############Output#######################
	  
	  IDMapping_output(mappingOutput,outputFileName,unMapF,dataType,inputGene,sourceIdType,targetIdType=standardID)
	 	r <- list(mapped=inputGene,unmapped=unMapF)
		return(r)
}


.processSourceIDMapPhosphosite <- function(hostName,organism,largeIdList,inputGene,standardID,dataType,idType,collapseMethod,methodType){
	
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
	if(is.character(mapR) && length(mapR)==1 && length(grep("ERROR:",mapR))>0){
		return(mapR)
	}
	
	mapF <- mapR$mapF
	unMapF <- mapR$unMapF
	
	colnames(mapF) <- c("phosphositeSeq",idType)
	
	idType_genelevel <- unlist(lapply(strsplit(mapF[,2],"_"),.combineG))
	idTypeGMap <- data.frame(phosi=mapF[,2],gene=idType_genelevel,glink="",stringsAsFactors=F)
	colnames(idTypeGMap)[1] <- idType
	
	#####Hard code#######
	if(length(grep("Uniprot",idType))>0){
		geneType <- "uniprot_swissprot"
		outLink <- "http://www.uniprot.org/uniprot/"
	}
	
	if(length(grep("Ensembl",idType))>0){
		geneType <- "ensembl_peptide_id"
		outLink <- paste("http://www.ensembl.org/",organism,"/Gene/Summary?db=core;t=",sep="")
	}
	
	if(length(grep("Refseq",idType))>0){
		geneType <- "refseq_peptide"
		outLink <- "https://www.ncbi.nlm.nih.gov/protein/"
	}
	idTypeGMap[,"glink"] <- paste(outLink,idTypeGMap[,2],sep="")
	
	######################
	
	########Get gene level information#########
	
	mapR <- IDMapping_gene(organism=organism,dataType="list",inputGene=unique(idTypeGMap[,"gene"]),sourceIdType=geneType,standardID="entrezgene",targetIdType="entrezgene",methodType=methodType,mappingOutput=FALSE,hostName=hostName)
	
	if(is.character(mapR) && length(mapR)==1 && length(grep("ERROR:",mapR))>0){
		return(mapR)
	}
	
	mapFG <- mapR$mapped
	mapFG <- mapFG[,-ncol(mapFG)]
	
	colnames(mapFG)[1] <- "gene"
	mergeA <- merge(x=idTypeGMap,y=mapFG,by="gene",all.x=TRUE)
	mergeA <- mergeA[,c(idType,"genesymbol","genename","glink")]
	mergeB <- merge(x=mapF,y=mergeA,by=idType,all.x=TRUE)
	
	
	if(dataType=="list"){
		inputGene <- mergeB[,c(idType,"genesymbol","genename",standardID,"glink")]
		colnames(inputGene)[1] <- "userid"
	}
	
	if(dataType=="rnk"){
		colnames(inputGene) <- c(idType,"score")
		inputGene <- merge(x=mergeB,y=inputGene,by=idType,all.x=TRUE)
	  inputGene <- inputGene[,c(idType,"genesymbol","genename",standardID,"score","glink")]
	  colnames(inputGene)[1] <- "userid"
	}
	
	if(dataType=="gmt"){
		colnames(inputGene) <- c("geneset","link",idType)
	  inputGene <- merge(x=mergeB,y=inputGene,by=idType,all.x=TRUE)
	  inputGene <- as.matrix(inputGene)
	  inputGene <- inputGene[,c("geneset","link",idType,"genesymbol","genename",standardID,"glink")]
	  colnames(inputGene)[3] <- "userid"
	}
	

	re <- list(mapped=inputGene,unmapped=unMapF)
	return(re)
}

.combineG <- function(e){
	e <- e[-length(e)]
	e <- paste(e,collapse="_")
	return(e)
}
