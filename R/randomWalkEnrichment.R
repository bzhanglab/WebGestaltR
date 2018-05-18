randomWalkEnrichment <- function(organism, network, method, inputSeed, topRank, seedNum, sigMethod, fdrThr, topThr, projectDir, projectName, hostName) {
	fileName <- paste(projectName, network, method, sep=".")
	geneSetUrl <- file.path(hostName, "api", "geneset")
	response <- GET(geneSetUrl, query=list(organism=organism, database=network, standardid="entrezgene", filetype="net"))
	net <- as.matrix(fread(content(response), sep='\t', header=FALSE, showProgress=FALSE))
	netGraph <- graph.edgelist(net, directed=FALSE)
	netNode <- V(netGraph)$name

	gmtUrl <- modify_url(geneSetUrl, query=list(organism=organism, database="geneontology_Biological_Process", standardid="genesymbol", filetype="gmt"))
	goAnn <- readGmt(gmtUrl)

	cat("Start Random Walk...\n")

	seeds <- unlist(strsplit(inputSeed, ","))

	allNum <- length(seeds)
	write.table(seeds, file.path(projectDir, paste0(fileName, "_seeds.txt")), row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
	seeds <- intersect(netNode, seeds)
	write.table(seeds, file.path(projectDir, paste0(fileName, "_seedsInNetwork.txt")), row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

	if (method == "Network_Retrieval_Prioritization") {
		if (length(seeds) < seedNum) {
			seedNum <- length(seeds)
		}
	}

	pt1 <- .netwalker(seeds, netGraph, r=0.5)
	gS <- data.frame(name=netNode, score=pt1, per=1, stringsAsFactors=F)
	if (method == "Network_Expansion") {
		gS <- gS[gS[, "name"] %in% setdiff(gS[, "name"], seeds), ]
		gS <- gS[order(-gS[, "score"]), ]
		candidate <- gS[1:topRank, ]
		allN <- c(seeds, candidate[, "name"])
	} else {
		gS <- gS[gS[, "name"] %in% seeds, ]
		gS <- gS[order(-gS[, "score"]), ]
		highSeeds <- gS[1:seedNum, "name"]
		allN <- seeds
		candidate <- gS
	}

	subNet <- net[net[, 1] %in% allN & net[, 2] %in% allN, ]
	allN <- union(subNet[, 1], subNet[, 2])

	if (length(allN) != 0){
		termInfo <- .enrichmentFunction(organism, netNode, allN, goAnn, seeds, sigMethod, fdrThr, topThr, hostName)
	} else {
		termInfo <- NULL
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
	write.table(subNet, file.path(projectDir, paste0(fileName, "_randomWalkNetwork.txt")), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
	write.table(candidate, file.path(projectDir, paste0(fileName, "_candidate.txt")), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
	write.table(termInfo, file.path(projectDir, paste0(fileName, "_enrichedResult.txt")), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}


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


.enrichmentFunction <- function(organism, reference, interest, goAnn, seeds, sigMethod, fdrThr, topThr, hostName) {
	goAnn <- goAnn[, c(3, 1)]
	annRef <- goAnn[goAnn[, 1] %in% reference, ]
	annInterest <- goAnn[goAnn[, 1] %in% interest, ]

	allRefNum <- length(unique(annRef[, 1]))
	allInterestNum <- length(unique(annInterest[, 1]))

	refTermCount <- tapply(annRef[, 1], annRef[, 2], length)

	geneSetUrl <- file.path(hostName, "api", "geneset")
	response <- POST(geneSetUrl, body=list(organism=organism, database="geneontology_Biological_Process",
									filetype="des", ids=unique(annRef[, 2])), encode="json")
	refTermName <- fread(content(response), sep="\t", col.names=c("id", "name"), header=FALSE, data.table=FALSE)
	rownames(refTermName) <- refTermName[, 1]
	refTermName <- refTermName[names(refTermCount), ]

	refTermCount <- data.frame(goId=names(refTermCount), name=refTermName[, 2], refNum=refTermCount, stringsAsFactors=FALSE)
	refTermCount <- refTermCount[order(refTermCount[, 1]), ]
	interestTermCount <- tapply(annInterest[, 1], annInterest[, 2], length)

	interestTermGene <- tapply(annInterest[, 1], annInterest[, 2], .getGenes, seeds)


	interestTermCount <- data.frame(goId=names(interestTermCount), interestNum=interestTermCount, interestGene=interestTermGene, stringsAsFactors=FALSE)
	interestTermCount <- interestTermCount[order(interestTermCount[, 1]), ]

	refInterestTermCount <- refTermCount

	refInterestTermCount$interestNum <- array(0, dim=c(length(refInterestTermCount$goId), 1))
	refInterestTermCount[refInterestTermCount$goId %in% interestTermCount[, 1], 4] <- interestTermCount$interestNum

	refInterestTermCount$interestGene <= array("", dim=c(length(refInterestTermCount$goId), 1))
	refInterestTermCount[refInterestTermCount$goId %in% interestTermCount[, 1], 5] <- interestTermCount$interestGene

	refInterestTermCount$expected <- array(0, dim=c(length(refInterestTermCount$goId), 1))
	refInterestTermCount$expected <- (allInterestNum / allRefNum) * refInterestTermCount[, 3]

	refInterestTermCount$ratio <- array(0, dim=c(length(refInterestTermCount$goId), 1))
	refInterestTermCount$ratio <- as.numeric(refInterestTermCount[, 4]) / as.numeric(refInterestTermCount[, 6])

	pv <- 1 - phyper(refInterestTermCount[, 4] - 1, allInterestNum, allRefNum - allInterestNum, refInterestTermCount[, 3], lower.tail=TRUE, log.p=FALSE)

	refInterestTermCount$pvalue <- pv
	adp <- p.adjust(pv, method="BH")
	refInterestTermCount$fdr <- adp
	refInterestTermCount <- refInterestTermCount[order(refInterestTermCount[, 9]), ]

	if(sigMethod=="fdr"){
		refInterestTermCountSig <- refInterestTermCount[refInterestTermCount$fdr < fdrThr, ]
	}else{
		refInterestTermCountSig <- refInterestTermCount[1:topThr, ]
	}
	refInterestTermCountSig <- refInterestTermCountSig[order(refInterestTermCountSig[, 9], refInterestTermCountSig[, 8]), ]

	return(refInterestTermCountSig)
}


.getGenes <- function(genelist, seeds) {
	genelist <- data.frame(gene=genelist, id=0, stringsAsFactors=F)
	genelist[genelist[, 1] %in% seeds, 2] <- 1
	genelist <- paste(genelist[, 1], genelist[, 2], sep="|")
	genelist <- paste(genelist, collapse=";")
	return(genelist)
}

