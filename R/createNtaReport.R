createNtaReport <- function(networkName, method, sigMethod, fdrThr, topThr, highlightOption, outputDirectory, projectDir, projectName, hostName) {
	namePrefix <- paste(projectName, networkName, method, sep=".")
	seedsFn <- file.path(projectDir, paste0(namePrefix, "_seedsInSubnetwork.txt"))
	networkFn <- file.path(projectDir, paste0(namePrefix, "_randomWalkNetwork.txt"))
	candidateFn <- file.path(projectDir, paste0(namePrefix, "_candidate.txt"))
	enrichFn <- file.path(projectDir, paste0(namePrefix, "_enrichedResult.txt"))
	summaryFn <- file.path(projectDir, paste0(namePrefix, "_resultSummary.txt"))
	jsonFn <- file.path(projectDir, paste0(namePrefix, ".json"))

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
	candidates <- fread(candidateFn, header=FALSE, col.names=c("candidate", "score", "label"))
	candidates$score <- sprintf("%2.2E", candidates$score)
	candidates <- unname(rowSplit(candidates))

	enrichment <- fread(enrichFn, header=FALSE, col.names=c("goId", "goName", "c", "o", "geneInfo", "expect", "ratio", "rawP", "adjP"))

	summary <- readLines(summaryFn)

	network <- fread(networkFn, header=FALSE, col.names=c("source", "target"))
	allNodes <- unique(c(network, recursive=TRUE))

	## Prepare JSON data for cytoscape of gene network, add edges data first
	cyJson <- lapply(split(network, seq(nrow(network))), function(x) list(data=x))
	## Highlight (larger node, colored red in table for expansion)
	## for expansion, option to decide seeds or neighbors
	## for retrieval, highlight seeds
	if (method == "Network_Expansion") {
		seeds <- seeds
		highlight <- sapply(allNodes, function(x) x %in% seeds)
		if (highlightOption != "Seeds") {
			highlight <- !highlight
		}
	} else {
		seeds <- scan(highSeedsFn, "character")
		highlight <- sapply(allNodes, function(x) x %in% seeds)
	}
	## add node data
	cyJson <- c(cyJson, mapply(function(x, y) list(data=list(id=x, highlight=y)), allNodes, highlight, SIMPLIFY=FALSE))
	cyJson <- toJSON(unname(cyJson))

	## Prepare GO enrichment table data
	geneStr <- vector("list", nrow(enrichment))
	highlightedGenes <- vector("list", nrow(enrichment))
	for (i in 1:nrow(enrichment)) {
		geneInfoList <- unlist(strsplit(enrichment[[i, "geneInfo"]], ";", fixed=TRUE))
		geneData <- vector("list", length(geneInfoList))
		for (j in 1:length(geneInfoList)) {
			splitInfo <- unlist(strsplit(geneInfoList[[j]], "|", fixed=TRUE))
			geneName <- splitInfo[1]
			label <- splitInfo[2]
			if ((method == "Network_Expansion" && highlightOption == "Seeds" && label != 0) ||
					(method == "Network_Expansion" && highlightOption != "Seeds" && label == 0)) {
				highlight <- TRUE
			} else {
				highlight <- FALSE
			}
			geneData[[j]] <- list(geneName=geneName, highlight=highlight)
		}
			highlightedGenes[[i]] <- sapply(geneData, function(x) { if (x$highlight) { return(x$geneName) }})
			highlightedGenes[[i]] <- paste(highlightedGenes[[i]][!sapply(highlightedGenes[[i]], is.null)], collapse=",")
		geneStr[[i]] <- paste(sapply(geneData, function(x) x[["geneName"]]), collapse=",")
	}
	## rowSplit convert data frame to list of lists. But nested loops are not handled by the function
	enrichment$geneListStr <- geneStr
	enrichment$highlightedGenes <- highlightedGenes
	enrichment <- unname(rowSplit(enrichment))

	version <- packageVersion("WebGestaltR")
	netName <- paste(networkName, "net", sep="_")
	dagName <- paste(networkName, "dag", sep="_")

	data <- list(networkName=networkName, summary=summary, networkJson=cyJson, dagJson=dagJson,
				method=method, seeds=seeds, highlightIsSeeds=highlightOption=="Seeds",
				methodIsNetworkExpansion=method=="Network_Expansion",
				candidates=candidates, sigMethodIsFdr=sigMethod=="fdr", fdrThr=fdrThr, topThr=topThr,
				enrichment=enrichment, hostName=hostName,
				dagName=dagName, netName=netName,
				zipPath=paste0(projectName, ".zip"),
				toolboxNet=list(name=netName, nodes=allNodes),
				toolboxDag=list(name=dagName, nodes=goDataList)
				)
	partials <- list(toolbox=readLines(system.file("templates/toolbox.mustache", package="WebGestaltR")))
	template <- readLines(system.file("templates/networkContent.mustache", package="WebGestaltR"))
	content <- whisker.render(template, data, partials=partials)

	data <- list(hostName=hostName,
				networkName=networkName,
				version=version,
				networkContent=content
				)
	header <- readLines(system.file("templates/header.mustache", package="WebGestaltR"))
	footer <- readLines(system.file("templates/footer.mustache", package="WebGestaltR"))
	partials <- list(header=header, footer=footer)
	template <- readLines(system.file("templates/ntaTemplate.mustache", package="WebGestaltR"))
	outFn <- file.path(projectDir, paste0("Report_", projectName, ".html"))
	cat(whisker.render(template, data, partials=partials), file=outFn)
}
