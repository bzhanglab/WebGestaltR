WebGestaltR <- function(enrichMethod="ORA", organism="hsapiens", enrichDatabase="geneontology_Biological_Process", enrichDatabaseFile=NULL, enrichDatabaseType=NULL, enrichDatabaseDescriptionFile=NULL, interestGeneFile=NULL, interestGene=NULL, interestGeneType=NULL, collapseMethod="mean", referenceGeneFile=NULL, referenceGene=NULL, referenceGeneType=NULL, referenceSet=NULL, minNum=10, maxNum=500, sigMethod="fdr", fdrMethod="BH", fdrThr=0.05, topThr=10, reportNum=20, perNum=1000, isOutput=TRUE, outputDirectory=getwd(), projectName=NULL, dagColor="continuous", setCoverNum=10, networkConstructionMethod=NULL, neighborNum=10, highlightType="Seeds", highlightSeedNum=10, nThreads=1, hostName="http://www.webgestalt.org/",  ...) {
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
	if ('lNum' %in% names(extraArgs)) {
		warning("Parameter lNum is obsolete.\n")
	}
	if ('dNum' %in% names(extraArgs)) {
		cat("WARNING: Parameter dNum is deprecated and changed to reportNum.\n")
	}

	## TODO: add para test for NTA
	errorTest <- parameterErrorMessage(enrichMethod=enrichMethod, organism=organism, collapseMethod=collapseMethod, minNum=minNum, maxNum=maxNum, fdrMethod=fdrMethod, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, reportNum=reportNum, perNum=perNum, isOutput=isOutput, outputDirectory=outputDirectory, dagColor=dagColor, hostName=hostName)
	if(!is.null(errorTest)){
		return(errorTest)
	}

	if(is.null(projectName)){
		projectName <- as.character(as.integer(Sys.time()))
	}
	projectName <- gsub('%', '_', projectName, fixed=TRUE) # use for GOSlim summary file name, need to escape %
	if (enrichMethod == "ORA") {
		enrichR <- WebGestaltROra(organism=organism, enrichDatabase=enrichDatabase, enrichDatabaseFile=enrichDatabaseFile, enrichDatabaseType=enrichDatabaseType, enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile,  interestGeneFile=interestGeneFile, interestGene=interestGene, interestGeneType=interestGeneType, collapseMethod=collapseMethod, referenceGeneFile=referenceGeneFile, referenceGene=referenceGene, referenceGeneType=referenceGeneType, referenceSet=referenceSet, minNum=minNum, maxNum=maxNum, fdrMethod=fdrMethod, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, reportNum=reportNum, setCoverNum=setCoverNum, isOutput=isOutput, outputDirectory=outputDirectory, projectName=projectName, dagColor=dagColor, nThreads=nThreads, hostName=hostName)
	} else if (enrichMethod == "GSEA") {
		enrichR <- WebGestaltRGsea(organism=organism, enrichDatabase=enrichDatabase, enrichDatabaseFile=enrichDatabaseFile, enrichDatabaseType=enrichDatabaseType, enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile,  interestGeneFile=interestGeneFile, interestGene=interestGene, interestGeneType=interestGeneType, collapseMethod=collapseMethod, minNum=minNum, maxNum=maxNum, fdrMethod=fdrMethod, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, reportNum=reportNum, setCoverNum=setCoverNum, perNum=perNum, isOutput=isOutput, outputDirectory=outputDirectory, projectName=projectName, dagColor=dagColor, nThreads=nThreads, hostName=hostName)
	} else if (enrichMethod == "NTA") {
		enrichR <- WebGestaltRNta(organism=organism, network=enrichDatabase, method=networkConstructionMethod, neighborNum=neighborNum, highlightSeedNum=highlightSeedNum, inputSeed=interestGene, inputSeedFile=interestGeneFile, interestGeneType=interestGeneType, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, outputDirectory=outputDirectory, projectName=projectName, highlightType=highlightType, hostName=hostName)
	}

	return(enrichR)
}
