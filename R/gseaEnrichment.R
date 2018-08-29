gseaEnrichment <- function (hostName, outputDirectory, projectName, geneRankList, geneSet, collapseMethod="mean", minNum=10, maxNum=500, sigMethod="fdr", fdrThr=0.05, topThr=10, perNum=1000, lNum=20, isOutput=TRUE) {
	projectFolder <- file.path(outputDirectory, paste("Project_", projectName, sep=""))
	if (!dir.exists(projectFolder)) {
		dir.create(projectFolder)
	}

	geneRankList[, 1] <- as.character(geneRankList[, 1])
	sortedScores <- geneRankList[order(geneRankList[, 2], decreasing=TRUE), 2]
	geneSet[, 3] <- as.character(geneSet[, 3])

	geneSetName <- as.data.frame(unique(geneSet[, c(1,2)]))
	colnames(geneSetName) <- c("geneset", "link")
	effectiveGeneSet <- geneSet[geneSet[, 3] %in% geneRankList[, 1], , drop=FALSE]

	geneSetNum <- tapply(effectiveGeneSet[, 3], effectiveGeneSet[, 1],length)
	geneSetNum <- geneSetNum[geneSetNum>=minNum & geneSetNum<=maxNum]
	if (length(geneSetNum)==0) {
		error <- paste("ERROR: The number of annotated IDs for all functional categories are not from ", minNum," to ", maxNum," for the GSEA enrichment method.",sep="")
		cat(error)
		return(error)
	}

	# collapse rank list
	a <- tapply(geneRankList[, 2], geneRankList[, 1], collapseMethod, na.rm=TRUE)
	geneRankList <- data.frame(geneid=names(a), score=unname(a), stringsAsFactors=FALSE)

	gseaRnk <- file.path(projectFolder, paste("Project_", projectName, "_GSEA.rnk", sep=""))
	write.table(geneRankList, file=gseaRnk, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

	outputF <- file.path(projectFolder, paste0("Project_", projectName, "_GSEA/"))
	relativeF <- file.path(".", paste0("Project_", projectName, "_GSEA"))
	if (!dir.exists(outputF)) {
		dir.create(outputF)
	}

	inputDf <- prepareInputMatrixGsea(geneRankList, effectiveGeneSet)

	gseaRes <- swGsea(inputDf, thresh_type="val", perms=perNum,
										min_set_size=minNum, max_set_size=maxNum,
										nThreads=8, rng_seed=as.integer(format(Sys.time(), "%H%M%S"))
										)
	enrichRes <- gseaRes$Enrichment_Results[c("ES", "NES", "p_val", "fdr")]
	names(enrichRes)[names(enrichRes) == 'p_val'] <- 'PValue'
	names(enrichRes)[names(enrichRes) == 'fdr'] <- 'FDR'
	# TODO: handle errors

	if (sigMethod == "fdr") {
		sig <- enrichRes[enrichRes$FDR < fdrThr, ]
		insig <- enrichRes[enrichRes$FDR >= fdrThr, ]
	} else if (sigMethod == "top") {
		enrichRes <- enrichRes[order(enrichRes$FDR, enrichRes$PValue), ]
		if (nrow(enrichRes) > topThr) {
			sig <- enrichRes[1: topThr, ]
			insig <- enrichRes[(topThr+1):nrow(enrichRes), ]
		} else {
			sig <- enrichRes
			insig <- NULL
		}
	}
	numSig = nrow(sig)
	if (numSig == 0) {
		cat("No significant set is identified based on FDR ", fdrThr, "!", sep="")
		return()
	}

	sig$geneset <- rownames(sig)
	if (!is.null(insig)) {
		insig$geneset <- rownames(insig)
		insig$leadingEdgeNum <- unname(sapply(insig$geneset, function(geneSetId) {
			rsum <- gseaRes$Running_Sums[, geneSetId]
			maxIndex <- match(max(rsum), rsum)
			return(sum(gseaRes$Items_in_Set[[geneSetId]]$rank <= maxIndex))
		}))
	}
	sig <- merge(sig, geneSetName, by="geneset")
	sig$Size <- unname(sapply(sig$geneset, function(x) nrow(gseaRes$Items_in_Set[[x]])))
	sig$plotPath <- unname(sapply(sig$geneset, function(x) file.path(relativeF, paste0(x, ".png"))))

	leadingGeneNum <- vector("integer", numSig)
	leadingGenes <- vector("character", numSig)
	for (i in 1:numSig) {
		geneSetId <- sig[i, "geneset"]
		rsum <- gseaRes$Running_Sums[, geneSetId]
		maxIndex <- match(max(rsum), rsum)
		indexes <- gseaRes$Items_in_Set[[geneSetId]]$rank <= maxIndex
		leadingGeneNum[[i]] <- sum(indexes)
		leadingGenes[[i]] <- paste(rownames(gseaRes$Items_in_Set[[geneSetId]])[indexes], collapse=";")
		genes <- gseaRes$Items_in_Set[[geneSetId]]

		# Plot GSEA-like enrichment plot
		png(file.path(outputF, paste0(geneSetId, ".png")), bg="transparent")
		plot.new()
		par(fig=c(0, 1, 0.5, 1), mar=c(0, 5, 2, 3), new=T)
		plot(1:length(gseaRes$Running_Sums[, geneSetId]), gseaRes$Running_Sums[, geneSetId],
			type="l", main=paste0("Enrichment plot: ", geneSetId),
			xlab="", ylab="Enrichment Score", xaxt='n')
		abline(v=maxIndex, lty=3)
		par(fig=c(0, 1, 0.35, 0.5), mar=c(0, 5, 0, 3), new=T)
		plot(genes$rank, rep(1, nrow(genes)), type="h",
			xlim=c(1, length(sortedScores)), ylim=c(0, 1), axes=F, ann=F)
		par(fig=c(0, 1, 0, 0.35), mar=c(4, 5, 0, 3), new=T)
		plot(1:length(sortedScores), sortedScores, type="h",
			ylab="Ranked list metric", xlab="Rank in Ordered Dataset")
		abline(v=maxIndex, lty=3)
		dev.off()
	}
	sig$leadingEdgeNum <- leadingGeneNum
	sig$leadingEdgeID <- leadingGenes

	return(list(enriched=sig, background=insig))
}
