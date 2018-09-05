WebGestaltRNta <- function(organism="hsapiens", network="network_PPI_BIOGRID", method="Network_Retrieval_Prioritization", inputSeed, inputSeedFile, edgeNum=10, seedNum=10, sigMethod="fdr", fdrThr=0.05, topThr=10, highlightOption="Seeds", outputDirectory=getwd(), projectName=NULL, hostName="http://www.webgestalt.org/") {

	if(is.null(projectName)){
		projectName <- as.character(as.integer(Sys.time()))
	}
	projectDir <- file.path(outputDirectory, paste0("Project_", projectName))
	dir.create(projectDir)

	inputGene <- formatCheck("list", inputGeneFile=inputSeedFile, inputGene=inputSeed)

	geneSetUrl <- file.path(hostName, "api", "geneset")
	response <- GET(geneSetUrl, query=list(organism=organism, database="geneontology_Biological_Process", standardid="entrezgene", filetype="dag"))
	dagInfo <- read_tsv(content(response), col_names=c("source", "target"), col_types="cc")

	## networks <- unlist(strsplit(network, ",", fixed=TRUE))
	## May need to bring back analysis of multiple networks
	fileName <- paste(projectName, network, method, sep=".")
	randomWalkEnrichment(organism=organism, network=network, method=method, seedNum=seedNum, inputSeed=inputGene,
						sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, projectDir=projectDir,
						topRank=edgeNum, projectName=projectName, hostName=hostName)

	enrichResFile <- file.path(projectDir, paste0(fileName, "_enrichedResult.txt"))

	goEnrichData <- sapply(strsplit(readLines(enrichResFile), "\t", fixed=TRUE), function(x) x[[1]])
	goTermList <- goEnrichData

	queue <- as.list(goTermList)
	## expand to include all linked nodes
	edges <- list()
	while (length(queue) > 0 ) {
		goTerm <- queue[[1]]
		if (length(queue) == 1) {
			queue <- list()
		} else {
			queue <- queue[2:length(queue)]
		}
		inEdges <- filter(dagInfo, target == goTerm)
		if (nrow(inEdges) > 0) {
			edges <- c(edges, lapply(split(inEdges, seq(nrow(inEdges))), function(x) list("data"=x)))
		}
		for (parentNode in inEdges$source) {
			if (!parentNode %in% goTermList) {
				queue <- c(queue, parentNode)
				goTermList <- c(goTermList, parentNode)
			}
		}
	}

	response <- POST(geneSetUrl, body=list(organism=organism, database="geneontology_Biological_Process",
		filetype="des", ids=goTermList), encode="json")
	goId2Term <- read_tsv(content(response), col_names=c("id", "name"), col_types="cc")

	jsonFile <- file.path(projectDir, paste0(fileName, ".json"));
	jsonData <- vector(mode="list", length=length(goTermList))
	inputEndIndex <- length(goEnrichData)
	for (i in 1:length(goTermList)) {
		goId <- goTermList[[i]]
		goName <- filter(goId2Term, id == goId)[[1, "name"]]
		dataSets <- i <= inputEndIndex
		jsonData[[i]] <- list(data=list(id=goId, name=goName, datasets=dataSets))
	}
	jsonData <- unname(c(jsonData, edges))

	cat(toJSON(jsonData), "\n", sep="", file=jsonFile)

	createNtaReport(networkName=network, method=method, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr,
					highlightOption=highlightOption, outputDirectory=outputDirectory, projectDir=projectDir,
					projectName=projectName, hostName=hostName)

	cwd <- getwd()
	setwd(projectDir)
	zip(paste0(projectName, ".zip"), ".", flags="-rq")
	setwd(cwd)

	cat("Results can be found in the ", projectDir, "!\n", sep="")

}
