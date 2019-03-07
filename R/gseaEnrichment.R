#' @importFrom dplyr select distinct filter arrange mutate left_join %>%
#' @importFrom readr write_tsv
gseaEnrichment <- function (hostName, outputDirectory, projectName, geneRankList, geneSet, geneSetDes=NULL, collapseMethod="mean", minNum=10, maxNum=500, sigMethod="fdr", fdrThr=0.05, topThr=10, perNum=1000, isOutput=TRUE, nThreads=1) {
	projectFolder <- file.path(outputDirectory, paste("Project_", projectName, sep=""))
	if (!dir.exists(projectFolder)) {
		dir.create(projectFolder)
	}

	colnames(geneRankList) <- c("gene", "score")
	sortedScores <- sort(geneRankList$score, decreasing=TRUE)

	geneSetName <- geneSet %>% select(.data$geneSet, link=.data$description) %>% distinct()
	effectiveGeneSet <- geneSet %>% filter(.data$gene %in% geneRankList$gene)

	geneSetNum <- tapply(effectiveGeneSet$gene, effectiveGeneSet$geneSet, length)
	geneSetNum <- geneSetNum[geneSetNum>=minNum & geneSetNum<=maxNum]
	if (length(geneSetNum)==0) {
		error <- paste("ERROR: The number of annotated IDs for all functional categories are not from ", minNum," to ", maxNum, " for the GSEA enrichment method.\n", sep="")
		cat(error)
		return(error)
	}

	# collapse rank list
	a <- tapply(geneRankList$score, geneRankList$gene, collapseMethod, na.rm=TRUE)
	geneRankList <- data.frame(gene=names(a), score=unname(a), stringsAsFactors=FALSE)

	gseaRnk <- file.path(projectFolder, paste("Project_", projectName, "_GSEA.rnk", sep=""))
	write_tsv(geneRankList, gseaRnk, col_names=FALSE)

	outputF <- file.path(projectFolder, paste0("Project_", projectName, "_GSEA/"))
	relativeF <- file.path(".", paste0("Project_", projectName, "_GSEA"))
	if (!dir.exists(outputF) && isOutput) {
		dir.create(outputF)
	}

	inputDf <- prepareInputMatrixGsea(geneRankList, effectiveGeneSet)

	gseaRes <- swGsea(inputDf, thresh_type="val", perms=perNum,
		min_set_size=minNum, max_set_size=maxNum,
		nThreads=nThreads, rng_seed=as.integer(format(Sys.time(), "%H%M%S"))
	)
	enrichRes <- gseaRes$Enrichment_Results %>%
		mutate(geneSet = rownames(gseaRes$Enrichment_Results)) %>%
		select(.data$geneSet, .data$ES, .data$NES, pValue=.data$p_val, FDR=.data$fdr)
	# TODO: handle errors

	if (sigMethod == "fdr") {
		sig <- filter(enrichRes, .data$FDR < fdrThr)
		insig <- filter(enrichRes, .data$FDR >= fdrThr)
	} else if (sigMethod == "top") {
		enrichRes <- arrange(enrichRes, .data$FDR, .data$pValue)
		tmpRes <- getTopGseaResults(enrichRes, topThr)
		sig <- tmpRes[[1]]
		insig <-tmpRes[[2]]
	}
	numSig = nrow(sig)
	if (numSig == 0) {
		error <- paste0("ERROR: No significant set is identified based on FDR ", fdrThr, "!\n")
		cat(error)
		return(error)
	}

	if (!is.null(insig)) {
		insig$leadingEdgeNum <- unname(sapply(insig$geneSet, function(geneSet) {
			rsum <- gseaRes$Running_Sums[, geneSet] # Running sum is a matrix of gene by gene set
			maxPeak <- max(rsum)
			minPeak <- min(rsum)
			if (abs(maxPeak) >= abs(minPeak)) {
				peakIndex <- match(max(rsum), rsum)
				leadingEdgeNum <- sum(gseaRes$Items_in_Set[[geneSet]]$rank <= peakIndex)
			} else {
				peakIndex <- match(min(rsum), rsum)
				leadingEdgeNum <- sum(gseaRes$Items_in_Set[[geneSet]]$rank >= peakIndex)
			}
			return(leadingEdgeNum)
		}))
	}
	sig <- sig %>% left_join(geneSetName, by="geneSet") %>%
		mutate(size = unname(sapply(geneSet, function(x) nrow(gseaRes$Items_in_Set[[x]])))) %>%
		mutate(plotPath = unname(sapply(geneSet, function(x) file.path(relativeF, paste0(sanitizeFileName(x), ".png")))))

	leadingGeneNum <- vector("integer", numSig)
	leadingGenes <- vector("character", numSig)
	for (i in 1:numSig) {
		geneSet <- sig[[i, "geneSet"]]
		es <- sig[[i, "ES"]]
		genes <- gseaRes$Items_in_Set[[geneSet]] # rowname is gene and one column called rank
		rsum <- gseaRes$Running_Sums[, geneSet]
		peakIndex <- match(ifelse(es > 0, max(rsum), min(rsum)), rsum)
		if (es > 0) {
			indexes <- genes$rank <= peakIndex
		} else {
			indexes <- genes$rank >= peakIndex
		}
		leadingGeneNum[[i]] <- sum(indexes)
		leadingGenes[[i]] <- paste(rownames(genes)[indexes], collapse=";")

		if (isOutput) {
			# Plot GSEA-like enrichment plot
			if (!is.null(geneSetDes)) {
				# same name of variable and column name, use quasiquotation !!
				title <- as.character((geneSetDes %>% filter(.data$geneSet == !!geneSet))[1, "description"])
			} else {
				title <- geneSet
			}
			wrappedTitle <- strwrap(paste0("Enrichment plot: ", title), 60)
			png(file.path(outputF, paste0(sanitizeFileName(geneSet), ".png")), bg="transparent", width=2000, height=2000)
			plot.new()
			par(fig=c(0, 1, 0.5, 1), mar=c(0, 6, 6 * length(wrappedTitle), 2), cex.axis=2.5, cex.main=5, cex.lab=3.2, lwd=2, new=TRUE)
			plot(1:length(gseaRes$Running_Sums[, geneSet]), gseaRes$Running_Sums[, geneSet],
				type="l", main=paste(wrappedTitle, collapse="\n"),
				xlab="", ylab="Enrichment Score", xaxt='n', lwd=3)
			abline(v=peakIndex, lty=3)
			par(fig=c(0, 1, 0.35, 0.5), mar=c(0, 6, 0, 2), new=TRUE)
			plot(genes$rank, rep(1, nrow(genes)), type="h",
				xlim=c(1, length(sortedScores)), ylim=c(0, 1), axes=FALSE, ann=FALSE)
			par(fig=c(0, 1, 0, 0.35), mar=c(6, 6, 0, 2), cex.axis=2.5, cex.lab=3.2, new=TRUE)
			plot(1:length(sortedScores), sortedScores, type="h",
				ylab="Ranked list metric", xlab="Rank in Ordered Dataset")
			abline(v=peakIndex, lty=3)
			dev.off()
		}
	}
	sig$leadingEdgeNum <- leadingGeneNum
	sig$leadingEdgeId <- leadingGenes

	return(list(enriched=sig, background=insig))
}
