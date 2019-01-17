#' @importFrom httr GET POST content
#' @importFrom readr read_tsv cols
#' @importFrom rjson toJSON
WebGestaltRNta <- function(organism="hsapiens", network="network_PPI_BIOGRID", method="Network_Retrieval_Prioritization", inputSeed, inputSeedFile, interestGeneType="genesymbol", neighborNum=10, highlightSeedNum=10, sigMethod="fdr", fdrThr=0.05, topThr=10, highlightType="Seeds", outputDirectory=getwd(), projectName=NULL, hostName="http://www.webgestalt.org/") {
	projectDir <- file.path(outputDirectory, paste0("Project_", projectName))
	dir.create(projectDir)

	inputGene <- formatCheck("list", inputGeneFile=inputSeedFile, inputGene=inputSeed)
	# only networks are in gene symbols
	# mapping always returns gene symbol, could map to genesymbol but takes two requests
	inputGene <- idMappingGene(organism=organism, dataType="list", inputGene=inputGene, sourceIdType=interestGeneType, targetIdType="entrezgene", mappingOutput=FALSE, hostName=hostName)
	inputGene <- inputGene$mapped$geneSymbol

	geneSetUrl <- file.path(hostName, "api", "geneset")
	response <- GET(geneSetUrl, query=list(organism=organism, database="geneontology_Biological_Process", standardId="entrezgene", fileType="dag"))
	dagInfo <- read_tsv(content(response), col_names=c("source", "target"), col_types="cc")

	## networks <- unlist(strsplit(network, ",", fixed=TRUE))
	## May need to bring back analysis of multiple networks
	fileName <- paste(projectName, network, method, sep=".")
	randomWalkEnrichment(organism=organism, network=network, method=method, highlightSeedNum=highlightSeedNum, inputSeed=inputGene,
						sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, projectDir=projectDir,
						topRank=neighborNum, projectName=projectName, hostName=hostName)

	enrichResFile <- file.path(projectDir, paste0(fileName, "_enrichedResult.txt"))

	goTermList <- read_tsv(enrichResFile, col_types=cols())$goId
	inputEndIndex <- length(goTermList)

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
		inEdges <- filter(dagInfo, .data$target == goTerm)
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
		fileType="des", ids=goTermList), encode="json")
	goId2Term <- read_tsv(content(response), col_names=c("id", "name"), col_types="cc")

	jsonFile <- file.path(projectDir, paste0(fileName, ".json"));
	jsonData <- vector(mode="list", length=length(goTermList))

	for (i in 1:length(goTermList)) {
		goId <- goTermList[[i]]
		goName <- filter(goId2Term, .data$id == goId)[[1, "name"]]
		dataSets <- i <= inputEndIndex
		jsonData[[i]] <- list(data=list(id=goId, name=goName, datasets=dataSets))
	}
	jsonData <- unname(c(jsonData, edges))

	cat(toJSON(jsonData), "\n", sep="", file=jsonFile)

	createNtaReport(networkName=network, method=method, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr,
					highlightType=highlightType, outputDirectory=outputDirectory, projectDir=projectDir,
					projectName=projectName, hostName=hostName)

	cwd <- getwd()
	setwd(projectDir)
	zip(paste0(projectName, ".zip"), ".", flags="-rq")
	setwd(cwd)

	cat("Results can be found in the ", projectDir, "!\n", sep="")

}
