WebGestaltR <- function(enrichMethod="ORA", organism="hsapiens", enrichDatabase="geneontology_Biological_Process", enrichDatabaseFile=NULL, enrichDatabaseType=NULL, enrichDatabaseDescriptionFile=NULL, interestGeneFile=NULL, interestGene=NULL, interestGeneType=NULL, collapseMethod="mean", referenceGeneFile=NULL, referenceGene=NULL, referenceGeneType=NULL, referenceSet=NULL, minNum=10, maxNum=500, fdrMethod="BH", sigMethod="fdr", fdrThr=0.05, topThr=10, dNum=20, perNum=1000, lNum=20, isOutput=TRUE, outputDirectory=getwd(), projectName=NULL, keepGseaFolder=FALSE, dagColor="continuous", hostName="http://www.webgestalt.org/",  ...) {

	extraArgs <- list(...)
	if ('keepGSEAFolder' %in% names(extraArgs)) {
		keepGseaFolder <- extraArgs$keepGSEAFolder
		cat("WARNING: Parameter keepGSEAFolder is deprecated and changed to keepGseaFolder!\n")
	}
	if ('is.output' %in% names(extraArgs)) {
		isOutput <- extraArgs$is.output
		cat("WARNING: Parameter is.output is deprecated and changed to isOutput!\n")
	}
	if ('methodType' %in% names(extraArgs)) {
		cat("WARNING: Parameter methodType is obsolete.\n")
	}

	errorTest <- parameterErrorMessage(enrichMethod=enrichMethod, organism=organism, collapseMethod=collapseMethod, minNum=minNum, maxNum=maxNum, fdrMethod=fdrMethod, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, dNum=dNum, perNum=perNum, lNum=lNum, isOutput=isOutput, outputDirectory=outputDirectory, keepGseaFolder=keepGseaFolder, dagColor=dagColor, hostName=hostName)
	if(!is.null(errorTest)){
		return(errorTest)
	}

	if(enrichMethod=="ORA"){
		enrichR <- WebGestaltROra(enrichMethod="ORA", organism=organism, enrichDatabase=enrichDatabase, enrichDatabaseFile=enrichDatabaseFile, enrichDatabaseType=enrichDatabaseType, enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile,  interestGeneFile=interestGeneFile, interestGene=interestGene, interestGeneType=interestGeneType, collapseMethod=collapseMethod, referenceGeneFile=referenceGeneFile, referenceGene=referenceGene, referenceGeneType=referenceGeneType, referenceSet=referenceSet, minNum=minNum, maxNum=maxNum, fdrMethod=fdrMethod, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, dNum=dNum, isOutput=isOutput, outputDirectory=outputDirectory, projectName=projectName, dagColor=dagColor, hostName=hostName)
	}

	if(enrichMethod=="GSEA"){
		enrichR <- WebGestaltRGsea(enrichMethod="GSEA", organism=organism, enrichDatabase=enrichDatabase, enrichDatabaseFile=enrichDatabaseFile, enrichDatabaseType=enrichDatabaseType, enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile,  interestGeneFile=interestGeneFile, interestGene=interestGene, interestGeneType=interestGeneType, collapseMethod=collapseMethod, minNum=minNum, maxNum=maxNum, fdrMethod=fdrMethod, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, dNum=dNum, perNum=perNum, lNum=lNum, isOutput=isOutput, outputDirectory=outputDirectory, projectName=projectName, keepGseaFolder=keepGseaFolder, dagColor=dagColor, hostName=hostName)
	}

	return(enrichR)
}
