WebGestaltR <- function(enrichMethod="ORA",organism="hsapiens",enrichDatabase="geneontology_Biological_Process",enrichDatabaseFile=NULL,enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL,interestGene=NULL,interestGeneType=NULL,collapseMethod="mean",referenceGeneFile=NULL,referenceGene=NULL,referenceGeneType=NULL,referenceSet=NULL,minNum=10,maxNum=500,fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,topThr=10,dNum=20,perNum=1000,lNum=20,is.output=TRUE,outputDirectory=getwd(),projectName=NULL,keepGSEAFolder=FALSE,methodType="R",dagColor="continuous",hostName="http://www.webgestalt.org/"){
		
		errorTest <- parameterErrorMessage(enrichMethod=enrichMethod,organism=organism,collapseMethod=collapseMethod,minNum=minNum,maxNum=maxNum,fdrMethod=fdrMethod,sigMethod=sigMethod,fdrThr=fdrThr,topThr=topThr,dNum=dNum,perNum=perNum,lNum=lNum,is.output=is.output,outputDirectory=outputDirectory,keepGSEAFolder=keepGSEAFolder,methodType=methodType,dagColor=dagColor,hostName=hostName)
		if(!is.null(errorTest)){
			return(errorTest)
		}
		
		if(enrichMethod=="ORA"){
			enrichR <- WebGestaltR_ORA(enrichMethod="ORA",organism=organism,enrichDatabase=enrichDatabase,enrichDatabaseFile=enrichDatabaseFile,enrichDatabaseType=enrichDatabaseType,enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile, interestGeneFile=interestGeneFile,interestGene=interestGene,interestGeneType=interestGeneType,collapseMethod=collapseMethod,referenceGeneFile=referenceGeneFile,referenceGene=referenceGene,referenceGeneType=referenceGeneType,referenceSet=referenceSet,minNum=minNum,maxNum=maxNum,fdrMethod=fdrMethod,sigMethod=sigMethod,fdrThr=fdrThr,topThr=topThr,dNum=dNum,is.output=is.output,outputDirectory=outputDirectory,projectName=projectName,methodType=methodType,dagColor=dagColor,hostName=hostName)
		}
		
		if(enrichMethod=="GSEA"){
			enrichR <- WebGestaltR_GSEA(enrichMethod="GSEA",organism=organism,enrichDatabase=enrichDatabase,enrichDatabaseFile=enrichDatabaseFile,enrichDatabaseType=enrichDatabaseType,enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile, interestGeneFile=interestGeneFile,interestGene=interestGene,interestGeneType=interestGeneType,collapseMethod=collapseMethod,minNum=minNum,maxNum=maxNum,fdrMethod=fdrMethod,sigMethod=sigMethod,fdrThr=fdrThr,topThr=topThr,dNum=dNum,perNum=perNum,lNum=lNum,is.output=is.output,outputDirectory=outputDirectory,projectName=projectName,keepGSEAFolder=keepGSEAFolder,methodType=methodType,dagColor=dagColor,hostName=hostName)
		}
		
		return(enrichR)
		
}