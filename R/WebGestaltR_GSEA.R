WebGestaltR_GSEA <- function(enrichMethod="GSEA",organism="hsapiens",enrichDatabase="geneontology_Biological_Process",enrichDatabaseFile=NULL,enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL, interestGeneFile=NULL,interestGene=NULL,interestGeneType=NULL,collapseMethod="mean",minNum=10,maxNum=500,fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,topThr=10,dNum=20,perNum=1000,lNum=20,is.output=TRUE,outputDirectory=getwd(),projectName=NULL,keepGSEAFolder=FALSE,methodType="R",dagColor="binary",hostName="http://www.webgestalt.org/"){
		
		
		if(is.null(projectName)){
			timeStamp <- gsub("\\.","_",as.character(as.numeric(Sys.time())))
    }else{
    	timeStamp <- projectName
    }
   
    projectDir <- file.path(outputDirectory,paste("Project_",timeStamp,sep=""))
    
		
		#########Web server will input "NULL" to the R package, thus, we need to change "NULL" to NULL########
		
		enrichDatabase <- testNull(enrichDatabase)
		enrichDatabaseFile <- testNull(enrichDatabaseFile)
		enrichDatabaseType <- testNull(enrichDatabaseType)
		enrichDatabaseDescriptionFile <- testNull(enrichDatabaseDescriptionFile)
		interestGeneFile <- testNull(interestGeneFile)
		interestGene <- testNull(interestGene)
		interestGeneType <- testNull(interestGeneType)
		
		
		################Check parameter################
		errorTest <- parameterErrorMessage(enrichMethod=enrichMethod,organism=organism,collapseMethod=collapseMethod,minNum=minNum,maxNum=maxNum,fdrMethod=fdrMethod,sigMethod=sigMethod,fdrThr=fdrThr,topThr=topThr,dNum=dNum,perNum=perNum,lNum=lNum,is.output=is.output,outputDirectory=outputDirectory,keepGSEAFolder=keepGSEAFolder,methodType=methodType,dagColor=dagColor,hostName=hostName)
		
		if(!is.null(errorTest)){
			return(errorTest)
		}
    
    ####################################
   	
		
		#############Check enriched database#############
		
		cat("Uploading the functional categories...\n")
		enrichD <- loadGeneSet(organism=organism,enrichDatabase=enrichDatabase,enrichDatabaseFile=enrichDatabaseFile,enrichDatabaseType=enrichDatabaseType,enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile,collapseMethod=collapseMethod,methodType=methodType,hostName=hostName)
		if(is.character(enrichD) && length(enrichD)==1 && length(grep("ERROR:",enrichD))>0){
   		return(enrichD)
   	}
   	
		geneSet <- enrichD$geneSet
   	geneSetDes <- enrichD$geneSetDes
   	geneSetDAG <- enrichD$geneSetDAG
   	geneSetNet <- enrichD$geneSetNet
   	database_standardId <- enrichD$standardId
		
		###########Check input interesting gene list###############
    cat("Uploading the ID list...\n")
    interestingGeneMap <- loadInterestGene(organism=organism,dataType="rnk",inputGeneFile=interestGeneFile,inputGene=interestGene,geneType=interestGeneType,collapseMethod=collapseMethod,methodType=methodType,hostName=hostName,geneSet=geneSet)
    
    if(is.character(interestingGeneMap) && length(interestingGeneMap)==1 && length(grep("ERROR:",interestingGeneMap))>0){
   		return(interestingGeneMap)
   	}
    
    if(organism=="others"){
    	interestGene_List <- unique(interestingGeneMap)
    }else{
    	interestStandardID <- interestingGeneMap$standardId
    	interestGene_List <- unique(interestingGeneMap$mapped[,c(interestStandardID,"score")])
    }
    	 
   ##########Create project folder##############
    if(is.output==TRUE){
     	dir.create(projectDir)
     
     	######Summary gene annotation based on the GOSlim###########
     
    	if(organism!="others"){
     		if(database_standardId=="entrezgene"){
     				cat("Summary the uploaded ID list by GO Slim data...\n")
     				goslim_output <- file.path(projectDir,paste("goslim_summary_",timeStamp,sep=""))
     				re <- GOSlimSummary(organism=organism,genelist=interestGene_List[,1],outputFile=goslim_output,outputType="png",hostName=hostName)
     				if(is.character(re) && length(re)==1 && length(grep("ERROR:",re))>0){
     					return(re)
     				}
     		}
     		
				write.table(interestingGeneMap$mapped,file=file.path(projectDir,paste("interestingID_Mappingtable_",timeStamp,".txt",sep="")),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
				write.table(interestingGeneMap$unmapped,file=file.path(projectDir,paste("interestingID_unmappedList_",timeStamp,".txt",sep="")),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
     	}else{
     		write.table(interestGene_List,file=file.path(projectDir,paste("interestList_",timeStamp,".txt",sep="")),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
     	}
    }
     
     
   #############Run enrichment analysis###################
		cat("Perform the enrichment analysis...\n")
		
		enrichedSig <- GSEAEnrichment(hostName,outputDirectory,timeStamp,interestGene_List,geneSet,minNum=minNum,maxNum=maxNum,sigMethod=sigMethod,fdrThr=fdrThr,topThr=topThr,perNum=perNum,lNum=lNum,is.output=is.output,keepGSEAFolder=keepGSEAFolder)
    
    if(is.character(enrichedSig) && length(enrichedSig)==1 && length(grep("ERROR:",enrichedSig))>0){
    	return(enrichedSig)
    }
    
    if(!is.null(enrichedSig)){
    	if(!is.null(geneSetDes)){ #######Add extra description information###########
    			colnames(geneSetDes) <- c("geneset","description")
    			enrichedSig <- merge(x=enrichedSig,y=geneSetDes,by="geneset",all.x=TRUE)
    			enrichedSig <- enrichedSig[,c(1,ncol(enrichedSig),2:(ncol(enrichedSig)-1))]
    			enrichedSig <- enrichedSig[order(enrichedSig[,"FDR"],enrichedSig[,"PValue"]),]
    	}
    	
    	if(organism!="others" && interestGeneType!=interestStandardID){
    			enrichedSig <- mapUserID(enrichedSig,"leadingEdgeID",interestingGeneMap)  ###mapUserID function is in the enrichmentResultProcess_component file
    	}
    	
    	if(is.output==TRUE){
    		write.table(enrichedSig,file=file.path(projectDir,paste("enrichment_results_",timeStamp,".txt",sep="")),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
    	}
    }
				
    if(is.output==TRUE){
    	
    	##############Create report##################
     	cat("Generate the final report...\n")
    	createReport(hostName=hostName,outputDirectory=outputDirectory,organism=organism,timeStamp=timeStamp,enrichMethod=enrichMethod,geneSet=geneSet,geneSetDes=geneSetDes,geneSetDAG=geneSetDAG,geneSetNet=geneSetNet,interestingGeneMap=interestingGeneMap,enrichedSig=enrichedSig,enrichDatabase=enrichDatabase,enrichDatabaseFile=enrichDatabaseFile,enrichDatabaseType=enrichDatabaseType,enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile,interestGeneFile=interestGeneFile,interestGene=interestGene,interestGeneType=interestGeneType,collapseMethod=collapseMethod,minNum=minNum,maxNum=maxNum,fdrMethod=fdrMethod,sigMethod=sigMethod,fdrThr=fdrThr,topThr=topThr,dNum=dNum,perNum=perNum,lNum=lNum,dagColor=dagColor)
    	
    	comm <- paste("tar -C ",projectDir," -zcvf ",projectDir,"/Project_",timeStamp,".tar.gz .",sep="")
    	system(comm,ignore.stderr=TRUE,ignore.stdout=TRUE)
    	
    	cat("Results can be found in the ",projectDir,"!",sep="")
    }
    return(enrichedSig)
}