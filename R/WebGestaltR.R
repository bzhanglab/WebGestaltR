WebGestaltR <- function(enrichMethod="ORA", organism="hsapiens", enrichDatabase="geneontology_Biological_Process", enrichDatabaseFile=NULL, enrichDatabaseType=NULL, enrichDatabaseDescriptionFile=NULL, interestGeneFile=NULL, interestGene=NULL, interestGeneType=NULL, collapseMethod="mean", referenceGeneFile=NULL, referenceGene=NULL, referenceGeneType=NULL, referenceSet=NULL, networkConstructionMethod=NULL, edgeNum=10, seedNum=10, minNum=10, maxNum=500, fdrMethod="BH", sigMethod="fdr", fdrThr=0.05, topThr=10, dNum=20, perNum=1000, lNum=20, isOutput=TRUE, outputDirectory=getwd(), projectName=NULL, dagColor="continuous", highlightOption="Seeds", hostName="http://www.webgestalt.org/",  ...) {

	extraArgs <- list(...)
	if ('keepGSEAFolder' %in% names(extraArgs) | 'keepGseaFolder' %in% names(extraArgs)) {
		cat("WARNING: Parameter keepGSEAFolder is obsolete.\n")
	}
	if ('is.output' %in% names(extraArgs)) {
		isOutput <- extraArgs$is.output
		cat("WARNING: Parameter is.output is deprecated and changed to isOutput!\n")
		warning("Column names of the result data frame are modified.")
	}
	if ('methodType' %in% names(extraArgs)) {
		cat("WARNING: Parameter methodType is obsolete.\n")
	}

	## TODO: add para test for NTA
	errorTest <- parameterErrorMessage(enrichMethod=enrichMethod, organism=organism, collapseMethod=collapseMethod, minNum=minNum, maxNum=maxNum, fdrMethod=fdrMethod, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, dNum=dNum, perNum=perNum, lNum=lNum, isOutput=isOutput, outputDirectory=outputDirectory, dagColor=dagColor, hostName=hostName)
	if(!is.null(errorTest)){
		return(errorTest)
	}

	if (enrichMethod == "ORA") {
		enrichR <- WebGestaltROra(organism=organism, enrichDatabase=enrichDatabase, enrichDatabaseFile=enrichDatabaseFile, enrichDatabaseType=enrichDatabaseType, enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile,  interestGeneFile=interestGeneFile, interestGene=interestGene, interestGeneType=interestGeneType, collapseMethod=collapseMethod, referenceGeneFile=referenceGeneFile, referenceGene=referenceGene, referenceGeneType=referenceGeneType, referenceSet=referenceSet, minNum=minNum, maxNum=maxNum, fdrMethod=fdrMethod, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, dNum=dNum, isOutput=isOutput, outputDirectory=outputDirectory, projectName=projectName, dagColor=dagColor, hostName=hostName)
	} else if (enrichMethod == "GSEA") {
		enrichR <- WebGestaltRGsea(organism=organism, enrichDatabase=enrichDatabase, enrichDatabaseFile=enrichDatabaseFile, enrichDatabaseType=enrichDatabaseType, enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile,  interestGeneFile=interestGeneFile, interestGene=interestGene, interestGeneType=interestGeneType, collapseMethod=collapseMethod, minNum=minNum, maxNum=maxNum, fdrMethod=fdrMethod, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, dNum=dNum, perNum=perNum, lNum=lNum, isOutput=isOutput, outputDirectory=outputDirectory, projectName=projectName, dagColor=dagColor, hostName=hostName)
	} else if (enrichMethod == "NTA") {
		enrichR <- WebGestaltRNta(organism=organism, network=enrichDatabase, method=networkConstructionMethod, edgeNum=edgeNum, seedNum=seedNum, inputSeed=interestGene, inputSeedFile=interestGeneFile, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, outputDirectory=outputDirectory, projectName=projectName, highlightOption=highlightOption, hostName=hostName)
	}

	return(enrichR)
}
