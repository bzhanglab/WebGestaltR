#' @importFrom httr POST content
#' @importFrom readr read_tsv cols
#' @importFrom jsonlite toJSON
WebGestaltRNta <- function(organism="hsapiens", network="network_PPI_BIOGRID", method="Network_Retrieval_Prioritization", inputSeed, inputSeedFile, interestGeneType="genesymbol", neighborNum=10, highlightSeedNum=10, sigMethod="fdr", fdrThr=0.05, topThr=10, highlightType="Seeds", outputDirectory=getwd(), projectName=NULL, cache=NULL, hostName="http://www.webgestalt.org/") {
	projectDir <- file.path(outputDirectory, paste0("Project_", projectName))
	dir.create(projectDir)

	if (length(network) > 1) {
		stop("NTA does not support multiple databases.")
	}
	inputGene <- formatCheck("list", inputGeneFile=inputSeedFile, inputGene=inputSeed)
	# only networks are in gene symbols
	# mapping always returns gene symbol, could map to genesymbol but takes two requests
	inputGene <- idMappingGene(organism=organism, dataType="list", inputGene=inputGene, sourceIdType=interestGeneType, targetIdType="entrezgene", mappingOutput=FALSE, hostName=hostName)
	inputGene <- inputGene$mapped$geneSymbol

	if (startsWith(hostName, "file://")) {
		dagInfo <- read_tsv(
			removeFileProtocol(file.path(hostName, "geneset", paste(organism, "geneontology_Biological_Process", "entrezgene.dag", sep="_"))),
			col_names=c("source", "target"), col_types="cc"
		)
	} else {
		geneSetUrl <- file.path(hostName, "api", "geneset")
		response <- cacheUrl(geneSetUrl, cache=cache, query=list(organism=organism, database="geneontology_Biological_Process", standardId="entrezgene", fileType="dag"))
		dagInfo <- read_tsv(content(response), col_names=c("source", "target"), col_types="cc")
	}

	## networks <- unlist(strsplit(network, ",", fixed=TRUE))
	## May need to bring back analysis of multiple networks
	fileName <- paste(projectName, network, method, sep=".")
	goEnrichRes <- randomWalkEnrichment(organism=organism, network=network, method=method, highlightSeedNum=highlightSeedNum, inputSeed=inputGene,
						sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr, projectDir=projectDir,
						topRank=neighborNum, projectName=projectName, cache=cache, hostName=hostName)
	if (is.null(goEnrichRes)) {
		return(NULL)
	}
	enrichResFile <- file.path(projectDir, paste0(fileName, "_enrichedResult.txt"))

	goTermList <- read_tsv(enrichResFile, col_types=cols())$goId
	inputEndIndex <- length(goTermList)

	dagTree <- expandDag(goTermList, dagInfo)
	goTermList <- dagTree$allNodes
	edges <- dagTree$edges
	rm(dagTree)

	if (startsWith(hostName, "file://")) {
		goId2Term <- read_tsv(
			removeFileProtocol(file.path(hostName, "geneset", paste(organism, "geneontology_Biological_Process", "entrezgene.des", sep="_"))),
			col_names=c("id", "name"), col_types="cc"
		)
	} else {
		response <- POST(geneSetUrl, body=list(organism=organism, database="geneontology_Biological_Process",
			fileType="des", ids=goTermList), encode="json")
		goId2Term <- read_tsv(content(response), col_names=c("id", "name"), col_types="cc")
	}

	jsonFile <- file.path(projectDir, paste0(fileName, ".json"));
	jsonData <- vector(mode="list", length=length(goTermList))

	for (i in 1:length(goTermList)) {
		goId <- goTermList[[i]]
		goName <- filter(goId2Term, .data$id == goId)[[1, "name"]]
		dataSets <- i <= inputEndIndex
		jsonData[[i]] <- list(data=list(id=goId, name=goName, datasets=dataSets))
	}
	jsonData <- unname(c(jsonData, edges))

	cat(toJSON(jsonData, auto_unbox=TRUE), "\n", sep="", file=jsonFile)

	createNtaReport(networkName=network, method=method, sigMethod=sigMethod, fdrThr=fdrThr, topThr=topThr,
					highlightType=highlightType, outputDirectory=outputDirectory, projectDir=projectDir,
					projectName=projectName, hostName=hostName)

	cwd <- getwd()
	setwd(projectDir)
	zip(paste0(projectName, ".zip"), ".", flags="-rq")
	setwd(cwd)

	cat("Results can be found in the ", projectDir, "!\n", sep="")
	return(goEnrichRes)
}
