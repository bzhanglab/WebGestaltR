#' Create HTML Report for NTA
#'
#' @importFrom readr read_tsv cols
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom whisker whisker.render rowSplit
#'
#' @keywords internal
#'
createNtaReport <- function(networkName, method, sigMethod, fdrThr, topThr, highlightType, outputDirectory, projectDir, projectName, hostName) {
	namePrefix <- paste(projectName, networkName, method, sep=".")
	seedsFn <- file.path(projectDir, paste0(namePrefix, "_seedsInSubnetwork.txt"))
	networkFn <- file.path(projectDir, paste0(namePrefix, "_randomWalkNetwork.txt"))
	candidateFn <- file.path(projectDir, paste0(namePrefix, "_candidate.txt"))
	enrichFn <- file.path(projectDir, paste0(namePrefix, "_enrichedResult.txt"))
	summaryFn <- file.path(projectDir, paste0(namePrefix, "_resultSummary.txt"))
	jsonFn <- file.path(projectDir, paste0(namePrefix, ".json"))

	if (startsWith(hostName, "file://")) {
		# change back hostName for web assets and browsers will cache it.
		hostName <- "https://www.webgestalt.org"
	}

	if (method == "Network_Retrieval_Prioritization") {
		highSeedsFn <- file.path(projectDir, paste0(namePrefix, "_highlightedSeeds.txt"))
	}
	## JSON for cytoscape of GO DAG
	dagJson <- readLines(jsonFn)
	dagData <- fromJSON(dagJson)
	## Get all GO term nodes, combine ID and Name
	goDataList <- lapply(dagData, function(x) paste(x[["data"]][["id"]], x[["data"]][["name"]]))
	goDataList <- as.character(goDataList[sapply(goDataList, length) > 0])

	seeds <- scan(seedsFn, "character")
	candidates <- read_tsv(candidateFn, col_names=c("candidate", "score"), col_types="cd-")
	candidates$score <- sprintf("%2.2E", candidates$score)

	enrichment <- read_tsv(enrichFn, col_names=c("goId", "goName", "c", "o", "geneInfo", "expect", "ratio", "rawP", "adjP"), skip=1, col_types=cols())

	summary <- readLines(summaryFn)

	network <- read_tsv(networkFn, col_names=c("source", "target"), col_types="cc")
	allNodes <- unique(c(network, recursive=TRUE))

	## Prepare JSON data for cytoscape of gene network, add edges data first
	cyJson <- lapply(split(network, seq(nrow(network))), function(x) list(data=as.list(x)))
	## Highlight (larger node, colored red in table for expansion)
	## for expansion, option to decide seeds or neighbors
	## for retrieval, highlight seeds
	if (method == "Network_Expansion") {
		seeds <- seeds
		highlight <- sapply(allNodes, function(x) x %in% seeds)
		if (highlightType != "Seeds") {
			highlight <- !highlight
		}
	} else {
		seeds <- scan(highSeedsFn, "character")
		highlight <- sapply(allNodes, function(x) x %in% seeds)
	}
	## add node data
	cyJson <- c(cyJson, mapply(function(x, y) list(data=list(id=x, highlight=y)), allNodes, highlight, SIMPLIFY=FALSE))
	cyJson <- toJSON(unname(cyJson), auto_unbox=TRUE)

	## Prepare GO enrichment table data
	#rowSplit still keeps each list element as data.frame
	enrichmentList <- lapply(unname(rowSplit(enrichment)), as.list)
	for (i in 1:nrow(enrichment)) {
		#decode geneInfo
		geneInfoList <- unlist(strsplit(enrichment[[i, "geneInfo"]], ";", fixed=TRUE))
		geneData <- vector("list", length(geneInfoList))
		for (j in 1:length(geneInfoList)) {
			splitInfo <- unlist(strsplit(geneInfoList[[j]], "|", fixed=TRUE))
			geneName <- splitInfo[1]
			label <- splitInfo[2]
			if ((method == "Network_Expansion" && highlightType == "Seeds" && label != 0) ||
					(method == "Network_Expansion" && highlightType != "Seeds" && label == 0)) {
				highlight <- TRUE
			} else {
				highlight <- FALSE
			}
			geneData[[j]] <- list(geneName=geneName, highlight=highlight)
		}
		enrichmentList[[i]]$geneInfo <- geneData
	}
	enrichmentJson <- toJSON(enrichmentList, auto_unbox=TRUE)

	version <- packageVersion("WebGestaltR")
	version <- paste(version[1, 1], version[1, 2], sep=".")
	netName <- paste(networkName, "net", sep="_")
	dagName <- paste(networkName, "dag", sep="_")

	data <- list(networkName=networkName, summary=summary, highlightType=highlightType,
				hostName=hostName,
				dagName=dagName, netName=netName,
				zipPath=paste0(projectName, ".zip"),
				toolboxDag=list(name=dagName, nodes=goDataList)
				)
	template <- readLines(system.file("templates/networkContent.mustache", package="WebGestaltR"))
	content <- whisker.render(template, data)

	data <- list(hostName=hostName,
				networkName=networkName,
				dagJson=dagJson,
				version=version,
				networkJson=cyJson,
				enrichmentJson=enrichmentJson,
				candidateJson=toJSON(candidates),
				seedJson=toJSON(seeds),
				method=method, sigMethod=sigMethod,
				threshold=ifelse(sigMethod=='fdr', fdrThr, topThr),
				networkContent=content
				)
	header <- readLines(system.file("templates/header.mustache", package="WebGestaltR"))
	footer <- readLines(system.file("templates/footer.mustache", package="WebGestaltR"))
	partials <- list(header=header, footer=footer)
	template <- readLines(system.file("templates/ntaTemplate.mustache", package="WebGestaltR"))
	outFn <- file.path(projectDir, paste0("Report_", projectName, ".html"))
	cat(whisker.render(template, data, partials=partials), file=outFn)

	file.remove(jsonFn)
}
