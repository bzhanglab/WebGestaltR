#' @importFrom httr GET content modify_url
#' @importFrom readr read_tsv
#' @importFrom dplyr filter arrange %>% desc
#' @importFrom igraph graph.edgelist V
randomWalkEnrichment <- function(organism, network, method, inputSeed, topRank, highlightSeedNum, sigMethod, fdrThr, topThr, projectDir, projectName, hostName) {
	fileName <- paste(projectName, network, method, sep=".")
	if (startsWith(hostName, "file://")) {
		net <- as.matrix(read_tsv(
			removeFileProtocol(file.path(hostName, "geneset", paste(organism, network, "entrezgene.net", sep="_"))),
			col_names=FALSE, col_types="cc"))
		goAnn <- readGmt(removeFileProtocol(file.path(hostName, "geneset", paste(organism, "geneontology_Biological_Process", "genesymbol.gmt", sep="_"))))
	} else {
		geneSetUrl <- file.path(hostName, "api", "geneset")
		# actually standard id is gene symbol for network
		response <- GET(geneSetUrl, query=list(organism=organism, database=network, standardId="entrezgene", fileType="net"))
		net <- as.matrix(read_tsv(content(response), col_names=FALSE, col_types="cc"))
		gmtUrl <- modify_url(geneSetUrl, query=list(organism=organism, database="geneontology_Biological_Process", standardId="genesymbol", fileType="gmt"))
		goAnn <- readGmt(gmtUrl)
	}
	netGraph <- graph.edgelist(net, directed=FALSE)
	netNode <- V(netGraph)$name


	cat("Start Random Walk...\n")

	#seeds <- unlist(strsplit(inputSeed, ","))
	seeds <- inputSeed

	allNum <- length(seeds)
	write(seeds, file.path(projectDir, paste0(fileName, "_seeds.txt")))
	seeds <- intersect(netNode, seeds)
	write(seeds, file.path(projectDir, paste0(fileName, "_seedsInNetwork.txt")))

	if (method == "Network_Retrieval_Prioritization") {
		if (length(seeds) < highlightSeedNum) {
			highlightSeedNum <- length(seeds)
		}
	}

	pt1 <- .netwalker(seeds, netGraph, r=0.5)
	gS <- data.frame(name=netNode, score=pt1, per=1, stringsAsFactors=F)
	if (method == "Network_Expansion") {
		gS <- gS %>% filter(.data$name %in% setdiff(.data$name, seeds)) %>%
			arrange(desc(.data$score))
		candidate <- gS[1:topRank, ]
		allN <- c(seeds, candidate$name)
	} else {
		gS <- gS %>% filter(.data$name %in% seeds) %>%
			arrange(desc(.data$score))
		highSeeds <- gS[1:highlightSeedNum, "name"]
		allN <- seeds
		candidate <- gS
	}

	subNet <- net[net[, 1] %in% allN & net[, 2] %in% allN, ]
	allN <- union(subNet[, 1], subNet[, 2])

	if (length(allN) != 0){
		termInfo <- .enrichmentFunction(organism, netNode, allN, goAnn, seeds, sigMethod, fdrThr, topThr, hostName)
	} else {
		stop("No sub-network is generated.")
	}

	cat("Output\n")
	if (method == "Network_Retrieval_Prioritization") {
		write(highSeeds, file.path(projectDir, paste0(fileName, "_highlightedSeeds.txt")))
	}
	overlapSeeds <- intersect(allN, seeds)
	x <- c(paste("Total number of genes in the selected network:", length(netNode), "(used for the enrichment analysis)"),
		paste("Total number of seeds:", allNum),
		paste("Total number of seeds in the selected network:", length(seeds))
		)
	if (method=="Network_Expansion") {
		x <- c(x,
			paste("Total number of genes in the expanded sub-network:", length(allN), "(used for the enrichment analysis)"),
			paste("Total number of seeds in the expanded sub-network:", length(overlapSeeds)),
			paste("We select top", topRank, "neighbors based on the probability of random walk method. All seeds and top ranking neighbors in the expanded sub-network can enrich to", nrow(termInfo), "GO BP categories.")
			)
	} else {
		x <- c(x,
			paste("Total number of seeds in the retrieved sub-network:", length(allN), "(used for the enrichment analysis)"),
			paste("All seeds in the retrieved sub-network can enrich to", nrow(termInfo), "GO BP categories.")
			)
	}

	write(x, file.path(projectDir, paste0(fileName, "_resultSummary.txt")))
	write(overlapSeeds, file.path(projectDir, paste0(fileName, "_seedsInSubnetwork.txt")))
	write_tsv(as.data.frame(subNet), file.path(projectDir, paste0(fileName, "_randomWalkNetwork.txt")), col_names=FALSE)
	write_tsv(candidate, file.path(projectDir, paste0(fileName, "_candidate.txt")), col_names=FALSE)
	if (!is.null(termInfo)) {
		write_tsv(termInfo, file.path(projectDir, paste0(fileName, "_enrichedResult.txt")))
	}
}

#' @importFrom igraph get.adjacency
.netwalker <- function(seed, network, r=0.5) {
	adjMatrix <- get.adjacency(network, sparse=FALSE)
	de <- apply(adjMatrix, 2, sum)
	w <- t(t(adjMatrix)/de)

	p0 <- array(0, dim=c(nrow(w), 1))
	rownames(p0) <- rownames(w)
	p0[seed, 1] <- 1/length(seed)

	pt <- p0
	pt1 <- (1-r)*(w%*%pt)+r*p0

	while (sum(abs(pt1-pt)) > 1e-6) {
		pt <- pt1
		pt1 <- (1-r)*(w%*%pt)+r*p0
	}

	return(pt1)
}

#' @importFrom dplyr select filter arrange left_join mutate %>%
#' @importFrom httr POST content
#' @importFrom readr read_tsv
#' @importFrom stats p.adjust phyper
.enrichmentFunction <- function(organism, reference, interest, goAnn, seeds, sigMethod, fdrThr, topThr, hostName) {
	goAnn <- select(goAnn, .data$gene, .data$geneSet)
	annRef <- filter(goAnn, .data$gene %in% reference)
	annInterest <- filter(goAnn, .data$gene %in% interest)

	allRefNum <- length(unique(annRef$gene))
	allInterestNum <- length(unique(annInterest$gene))

	refTermCount <- tapply(annRef$gene, annRef$geneSet, length)

	if (startsWith(hostName, "file://")) {
		refTermName <- read_tsv(
			removeFileProtocol(file.path(hostName, "geneset", paste(organism, "geneontology_Biological_Process", "entrezgene.des", sep="_"))),
			col_names=c("id", "name"), col_types="cc"
		) %>% filter(.data$id %in% names(refTermCount))
	} else {
		geneSetUrl <- file.path(hostName, "api", "geneset")
		response <- POST(geneSetUrl, body=list(organism=organism, database="geneontology_Biological_Process",
			fileType="des", ids=unique(annRef$geneSet)), encode="json")
		refTermName <- read_tsv(content(response), col_names=c("id", "name"), col_types="cc") %>%
			filter(.data$id %in% names(refTermCount))
	}

	refTermCount <- data.frame(goId=names(refTermCount), refNum=refTermCount, stringsAsFactors=FALSE) %>%
		left_join(refTermName, by=c("goId"="id")) %>%
		select(.data$goId, .data$name, .data$refNum) %>%
		arrange(.data$goId)
	interestTermCount <- tapply(annInterest$gene, annInterest$geneSet, length)

	interestTermGene <- tapply(annInterest$gene, annInterest$geneSet, .getGenes, seeds)


	interestTermCount <- data.frame(goId=names(interestTermCount), interestNum=interestTermCount, interestGene=interestTermGene, stringsAsFactors=FALSE) %>%
		arrange(.data$goId)

	refInterestTermCount <- refTermCount

	refInterestTermCount$interestNum <- vector("numeric", nrow(refInterestTermCount))
	refInterestTermCount[refInterestTermCount$goId %in% interestTermCount$goId, "interestNum"] <- interestTermCount$interestNum

	refInterestTermCount$interestGene <- vector("numeric", nrow(refInterestTermCount))
	refInterestTermCount[refInterestTermCount$goId %in% interestTermCount$goId, "interestGene"] <- interestTermCount$interestGene

	refInterestTermCount <- refInterestTermCount %>%
		mutate(expect = (allInterestNum / allRefNum) * .data$refNum,
			ratio = .data$interestNum / .data$expect,
			pValue = 1 - phyper(.data$interestNum - 1, allInterestNum, allRefNum - allInterestNum, .data$refNum, lower.tail=TRUE, log.p=FALSE),
			FDR = p.adjust(.data$pValue, method="BH")
		) %>%
		arrange(.data$FDR, .data$pValue)


	if(sigMethod=="fdr"){
		refInterestTermCountSig <- filter(refInterestTermCount, .data$FDR < fdrThr)
	}else{
		refInterestTermCountSig <- refInterestTermCount[1:topThr, ]
	}

	return(refInterestTermCountSig)
}


.getGenes <- function(genelist, seeds) {
	genelist <- data.frame(gene=genelist, id=0, stringsAsFactors=F)
	genelist[genelist[, 1] %in% seeds, 2] <- 1
	genelist <- paste(genelist[, 1], genelist[, 2], sep="|")
	genelist <- paste(genelist, collapse=";")
	return(genelist)
}

