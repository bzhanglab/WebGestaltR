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

	bpFileRe <- .processData(organism,hostName,geneList,"BiologicalProcess")
	if(.hasError(bpFileRe)){
		return(bpFileRe)
	}else{
		bpFile <- bpFileRe$goFile
		bpUnclassified <- bpFileRe$dataUnclassified
	}

	ccFileRe <- .processData(organism,hostName,geneList,"CellularComponent")
	if(.hasError(ccFileRe)){
		return(ccFileRe)
	}else{
		ccFile <- ccFileRe$goFile
		ccUnclassified <- ccFileRe$dataUnclassified
	}

	mfFileRe <- .processData(organism,hostName,geneList,"MolecularFunction")
	if(.hasError(mfFileRe)){
		return(mfFileRe)
	}else{
		mfFile <- mfFileRe$goFile
		mfUnclassified <- mfFileRe$dataUnclassified
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

	.plotData(geneList,bpFile,bpUnclassified,"Biological Process","red")
	.plotData(geneList,ccFile,ccUnclassified,"Cellular Component","blue")
	.plotData(geneList,mfFile,mfUnclassified,"Molecular Function","green")

	dev.off()
	return(NULL)
}

.processData <- function(organism,hostName,geneList,ontology){
	###file.path#######
	goFile <- file.path(hostName,"data","goslim",paste(organism,"_GOSlim_",ontology,".table",sep=""))

	goFile <- fread(input=goFile,header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses="character",data.table=FALSE)

	goFile <- goFile[goFile[,1] %in% geneList,,drop=FALSE]

	if(nrow(goFile)==0){
		error <- "ERROR: The ID type of the uploaded list is not entrez gene id. Please use IDMapping function to map other id types to entrez gene id."
		cat(error)
		return(error)
	}

	dataUnclassified <- setdiff(geneList,unique(goFile[,1]))

	goFile1 <- tapply(goFile[,1],goFile[,2],length)

	goFile1 <- goFile1[goFile1>0]

	goFile <- unique(goFile[,c(2,3)])

	goFile1 <- data.frame(goAcc=names(goFile1),geneNum=goFile1,stringsAsFactors=FALSE)

	colnames(goFile) <- c("goAcc","goName")

	goFile <- merge(x=goFile,y=goFile1,by="goAcc")

	goFile <- goFile[order(-goFile[,3]),]
	re <- list(goFile=goFile,dataUnclassified=dataUnclassified)
	return(re)
}

.plotData <- function(geneList,goFile,dataUnclassified,ontology,color){
	par(mar=c(20,5,2,2))
	c <- c(length(geneList),goFile[,3],length(dataUnclassified))
	names(c) <- c("all",goFile[,2],"unclassified")
	maxC <- max(c)
	xx <- barplot(c,main=paste("Bar chart of ",ontology," categories",sep=""),col=color,xaxt="n",xlab="",ylim=c(0,1.2*maxC),cex.axis=1.5,las=2,cex.main=1.7)
	text(x=xx+0.3,y=c+maxC*0.05,labels=c,pos=3,cex=1.5,col="black",srt=90)
	axis(side=1,at=xx,labels=names(c),las=2,cex.axis=1.5)
}


GOSlimSummary <- function(...) {
	cat("WARNING: Function GOSlimSummary is deprecated and changed to goSlimSummary!\n")
	return(goSlimSummary(...))
}
