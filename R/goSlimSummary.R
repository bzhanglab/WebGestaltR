goSlimSummary <- function(organism="hsapiens",geneList,outputFile,outputType="pdf",hostName="http://www.webgestalt.org"){
	organisms <- listOrganism(hostName=hostName)
	if(!(organism %in% organisms)){
		error <- paste("ERROR: ",organism," can not be supported.",sep="")
		cat(error)
		return(error)
	}

	if(!is.vector(geneList)){
		error <- "ERROR: Please upload a list of entrez gene ids.\n"
		cat(error)
		return(error)
	}else{
		geneList <- as.character(geneList)
	}

	outputTypeList <- c("pdf","png","bmp")
	if(!(outputType %in% outputTypeList)){
		error <- paste("ERROR: The output Type ",outputType," is invalid. Please select one from pdf, png, and bmp.",sep="")
		cat(error)
		return(error)
	}

	bpRes <- .processData(organism,hostName,geneList,"BiologicalProcess")
	if(.hasError(bpRes)){
		return(bpRes)
	}else{
		bpGoCnts<- bpRes$goTermCounts
		bpUnclassified <- bpRes$dataUnclassified
	}

	ccRes <- .processData(organism,hostName,geneList,"CellularComponent")
	if(.hasError(ccRes)){
		return(ccRes)
	}else{
		ccGoCnts <- ccRes$goTermCounts
		ccUnclassified <- ccRes$dataUnclassified
	}

	mfRes <- .processData(organism,hostName,geneList,"MolecularFunction")
	if(.hasError(mfRes)){
		return(mfRes)
	}else{
		mfGoCnts <- mfRes$goTermCounts
		mfUnclassified <- mfRes$dataUnclassified
	}

	if(outputType=="pdf"){
		pdf(paste(outputFile,".pdf",sep=""),height=8,width=16)
	}

	if(outputType=="png"){
		png(paste(outputFile,".png",sep=""),width=2800,height=1300,res=200)
	}

	if(outputType=="bmp"){
		bmp(paste(outputFile,".bmp",sep=""),width=2800,height=1300,res=200)
	}

	layout(matrix(1:3,1,3))

	.plotData(geneList,bpGoCnts,bpUnclassified,"Biological Process","red")
	.plotData(geneList,ccGoCnts,ccUnclassified,"Cellular Component","blue")
	.plotData(geneList,mfGoCnts,mfUnclassified,"Molecular Function","green")

	dev.off()
	return(NULL)
}

.processData <- function(organism,hostName,geneList,ontology){
	goUrl <- file.path(hostName, "api", "goslim")
	response <- POST(goUrl, body=list(organism=organism, ontology=ontology, entrezgenes=geneList), encode="json")
	if (response$status_code != 200) {
		return(webRequestError(response))
	}
	names <- c("entrezgene", "accession", "name")
	goSlimData <- content(response)[["goslim"]]
	goSlimData <- data.frame(matrix(unlist(lapply(goSlimData, FUN=function(x) { x[names] })), nrow=length(goSlimData), byrow=TRUE), stringsAsFactors=FALSE)
	colnames(goSlimData) <- names


	if(nrow(goSlimData)==0){
		error <- "ERROR: The ID type of the uploaded list is not entrez gene id. Please use IDMapping function to map other id types to entrez gene id."
		cat(error)
		return(error)
	}

	dataUnclassified <- setdiff(geneList,unique(goSlimData[["entrezgene"]]))

	goTermCount <- tapply(goSlimData[["entrezgene"]], goSlimData[["accession"]], length)
	goTermCount <- data.frame(accession=names(goTermCount),geneNum=goTermCount,stringsAsFactors=FALSE)
	uniqueGoTerms <- unique(goSlimData[, c("accession", "name")])

	uniqueGoCounts <- merge(x=uniqueGoTerms,y=goTermCount,by="accession")
	uniqueGoCounts <- uniqueGoCounts[order(-uniqueGoCounts[["geneNum"]]),]
	re <- list(goTermCounts=uniqueGoCounts,dataUnclassified=dataUnclassified)
	return(re)
}

.plotData <- function(geneList,goCounts,dataUnclassified,ontology,color){
	par(mar=c(20,5,2,2))
	c <- c(length(geneList),goCounts[,"geneNum"],length(dataUnclassified))
	names(c) <- c("all",goCounts[,"name"],"unclassified")
	maxC <- max(c)
	xx <- barplot(c,main=paste("Bar chart of ",ontology," categories",sep=""),col=color,xaxt="n",xlab="",ylim=c(0,1.2*maxC),cex.axis=1.5,las=2,cex.main=1.7)
	text(x=xx+0.3,y=c+maxC*0.05,labels=c,pos=3,cex=1.5,col="black",srt=90)
	axis(side=1,at=xx,labels=names(c),las=2,cex.axis=1.5)
}


GOSlimSummary <- function(...) {
	cat("WARNING: Function GOSlimSummary is deprecated and changed to goSlimSummary!\n")
	return(goSlimSummary(...))
}