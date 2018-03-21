WebGestaltRBatch <- function(interestGeneFolder=NULL,interestGeneType=NULL,enrichMethod="ORA",organism="hsapiens",enrichDatabase="geneontology_Biological_Process",enrichDatabaseFile=NULL,enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,collapseMethod="mean",referenceGeneFile=NULL,referenceGene=NULL,referenceGeneType=NULL,referenceSet=NULL,minNum=10,maxNum=500,fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,topThr=10,dNum=20,perNum=1000,lNum=20,isOutput=TRUE,outputDirectory=getwd(),keepGseaFolder=FALSE,methodType="R",hostName="http://www.webgestalt.org/",isParallel=FALSE,nThreads=3){

	###HARD CODE#####
	if(enrichMethod=="ORA"){
		interestGeneFiles <- list.files(interestGeneFolder,pattern="\\.txt",full.names=TRUE)
	}

	if(enrichMethod=="GSEA"){
		interestGeneFiles <- list.files(interestGeneFolder,pattern="\\.rnk",full.names=TRUE)
	}

	projectNames <- unlist(lapply(strsplit(basename(interestGeneFiles),split=".",fixed=TRUE),function(e){return(paste(e[-length(e)],collapse="."))}))

	resultList <- list()

	if(isParallel==TRUE){
		cl <- makeCluster(nThreads)
		registerDoParallel(cl)
		resultList <- foreach(i=1:length(interestGeneFiles), .packages="WebGestaltR") %dopar% {
			sig <- WebGestaltR(enrichMethod=enrichMethod,organism=organism,enrichDatabase=enrichDatabase,enrichDatabaseFile=enrichDatabaseFile,enrichDatabaseType=enrichDatabaseType,enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile,interestGeneFile=interestGeneFiles[i],interestGene=NULL,interestGeneType=interestGeneType,collapseMethod=collapseMethod,referenceGeneFile=referenceGeneFile,referenceGene=referenceGene,referenceGeneType=referenceGeneType,referenceSet=referenceSet,minNum=minNum,maxNum=maxNum,fdrMethod=fdrMethod,sigMethod=sigMethod,fdrThr=fdrThr,topThr=topThr,dNum=dNum,perNum=perNum,lNum=lNum,isOutput=isOutput,outputDirectory=outputDirectory,projectName=projectNames[i],keepGseaFolder=keepGseaFolder,methodType=methodType,hostName=hostName)
			re <- list(filename=interestGeneFiles[i],enrichResult=sig)
			return(re)
		}
		stopCluster(cl)
	}else{
		for(i in c(1:length(interestGeneFiles))){
			cat("Process file:",interestGeneFiles[i],"\n",sep="")
			sig <- WebGestaltR(enrichMethod=enrichMethod,organism=organism,enrichDatabase=enrichDatabase,enrichDatabaseFile=enrichDatabaseFile,enrichDatabaseType=enrichDatabaseType,enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile,interestGeneFile=interestGeneFiles[i],interestGene=NULL,interestGeneType=interestGeneType,collapseMethod=collapseMethod,referenceGeneFile=referenceGeneFile,referenceGene=referenceGene,referenceGeneType=referenceGeneType,referenceSet=referenceSet,minNum=minNum,maxNum=maxNum,fdrMethod=fdrMethod,sigMethod=sigMethod,fdrThr=fdrThr,topThr=topThr,dNum=dNum,perNum=perNum,lNum=lNum,isOutput=isOutput,outputDirectory=outputDirectory,projectName=projectNames[i],keepGseaFolder=keepGseaFolder,methodType=methodType,hostName=hostName)
			re <- list(filename=interestGeneFiles[i],enrichResult=sig)
			resultList[[i]] <- re
		}
	}
	return(resultList)
}

WebGestaltR_batch <- function(keepGSEAFolder=FALSE, is.output=TRUE, ...) {
	cat("WARNING: Function WebGestaltR_batch is deprecated and changed to WebGestaltRBatch!\n")
	return(WebGestaltRBatch(keepGseaFolder=keepGSEAFolder, isOutput=is.output, ...))
}
