GOSlimSummary <- function(organism="hsapiens",genelist,outputFile,outputType="pdf",hostName="http://www.webgestalt.org"){
	
	organisms <- listOrganism(hostName=hostName)
	if(length(which(organisms==organism))==0){
    	error <- paste("ERROR: ",organism," can not be supported.",sep="")
    	cat(error)
			return(error)
   }
   
   if(!is.vector(genelist)){
			error <- "ERROR: Please upload a list of entrez gene ids.\n"
			cat(error)
			return(error)
	 }else{
			genelist <- as.character(genelist)
	 }
   
   outputTypeList <- c("pdf","png","bmp")
   if(length(which(outputTypeList==outputType))==0){
   		error <- paste("ERROR: The output Type ",outputType," is invalid. Please select one from pdf, png, and bmp.",sep="")
			cat(error)
			return(error)
   }
   	
	bpfilere <- .processData(organism,hostName,genelist,"BiologicalProcess")
	if(is.character(bpfilere) && length(bpfilere)==1 && length(grep("ERROR:",bpfilere))>0){
		return(bpfilere)
	}else{
		bpfile <- bpfilere$goFile
		bp_unclassified <- bpfilere$data_unclassified
	}
	
	ccfilere <- .processData(organism,hostName,genelist,"CellularComponent")
	if(is.character(ccfilere) && length(ccfilere)==1 && length(grep("ERROR:",ccfilere))>0){
		return(ccfilere)
	}else{
		ccfile <- ccfilere$goFile
		cc_unclassified <- ccfilere$data_unclassified
	}
	
	mffilere <- .processData(organism,hostName,genelist,"MolecularFunction")
	if(is.character(mffilere) && length(mffilere)==1 && length(grep("ERROR:",mffilere))>0){
		return(mffilere)
	}else{
		mffile <- mffilere$goFile
		mf_unclassified <- mffilere$data_unclassified
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
	
	.plotData(genelist,bpfile,bp_unclassified,"Biological Process","red")
	.plotData(genelist,ccfile,cc_unclassified,"Cellular Component","blue")
	.plotData(genelist,mffile,mf_unclassified,"Molecular Function","green")
	
	dev.off()
	return(NULL)
}

.processData <- function(organism,hostName,genelist,ontology){
	
	###file.path#######
	goFile <- file.path(hostName,"data","goslim",paste(organism,"_GOSlim_",ontology,".table",sep=""))

	goFile <- fread(input=goFile,header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses="character",data.table=FALSE)
	
	goFile <- goFile[goFile[,1] %in% genelist,,drop=FALSE]

	if(nrow(goFile)==0){
		error <- "ERROR: The ID type of the uploaded list is not entrez gene id. Please use IDMapping function to map other id types to entrez gene id."
		cat(error)
		return(error)
	}
	
	data_unclassified <- setdiff(genelist,unique(goFile[,1]))
	
	goFile1 <- tapply(goFile[,1],goFile[,2],length)
	
	goFile1 <- goFile1[goFile1>0]
	
	goFile <- unique(goFile[,c(2,3)])
	
	goFile1 <- data.frame(goacc=names(goFile1),genenum=goFile1,stringsAsFactors=FALSE)
	
	colnames(goFile) <- c("goacc","goname")

	goFile <- merge(x=goFile,y=goFile1,by="goacc")
	
	goFile <- goFile[order(-goFile[,3]),]
	re <- list(goFile=goFile,data_unclassified=data_unclassified)
	return(re)
}

.plotData <- function(genelist,goFile,data_unclassified,ontology,color){
	par(mar=c(20,5,2,2))
	c <- c(length(genelist),goFile[,3],length(data_unclassified))
	names(c) <- c("all",goFile[,2],"unclassified")
	maxC <- max(c)
	xx <- barplot(c,main=paste("Bar chart of ",ontology," categories",sep=""),col=color,xaxt="n",xlab="",ylim=c(0,1.2*maxC),cex.axis=1.5,las=2,cex.main=1.7)
	text(x=xx+0.3,y=c+maxC*0.05,labels=c,pos=3,cex=1.5,col="black",srt=90)
	axis(side=1,at=xx,labels=names(c),las=2,cex.axis=1.5)
}