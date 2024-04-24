#' @importFrom dplyr select distinct filter arrange mutate left_join %>%
#' @importFrom readr write_tsv
multiGseaEnrichment <- function(hostName = NULL, outputDirectory = NULL, projectName = NULL, geneRankList_list = NULL, geneSet_list = NULL,
                                geneSetDes_list = NULL, collapseMethod = "mean", minNum = 10, maxNum = 500, sigMethod = "fdr", fdrThr = 0.05,
                                topThr = 10, perNum = 1000, p = 1, isOutput = TRUE, saveRawGseaResult = FALSE, plotFormat = "png", nThreads = 1,
                                listNames = NULL) {
    inputDf_list <- list()
    old_project_name <- projectName
    old_projectFolder <- file.path(outputDirectory, paste("Project_", old_project_name, sep = ""))
    if (!dir.exists(old_projectFolder)) {
        dir.create(old_projectFolder)
    }
    geneSetName_list <- list()
    for (i in seq_along(geneRankList_list)) {
        projectName <- paste0(old_project_name, "_", listNames[i])
        geneRankList <- geneRankList_list[[i]]
        geneSet <- geneSet_list[[i]]
        if (is.null(geneSetDes_list) || length(geneSetDes_list) == 0 || length(geneSetDes_list) < i) {
            geneSetDes <- NULL
        } else {
            geneSetDes <- geneSetDes_list[[i]]
        }
        projectFolder <- file.path(outputDirectory, paste("Project_", old_project_name, "/", projectName, sep = ""))
        if (!dir.exists(projectFolder)) {
            dir.create(projectFolder)
        }
        colnames(geneRankList) <- c("gene", "score")
        sortedScores <- sort(geneRankList$score, decreasing = TRUE)
        geneSetName <- geneSet %>%
            select(.data$geneSet, link = .data$description) %>%
            distinct(.keep_all = TRUE)
        geneSetName_list[[i]] <- geneSetName
        effectiveGeneSet <- geneSet %>% filter(.data$gene %in% geneRankList$gene)

        geneSetNum <- tapply(effectiveGeneSet$gene, effectiveGeneSet$geneSet, length)
        geneSetNum <- geneSetNum[geneSetNum >= minNum & geneSetNum <= maxNum]
        if (length(geneSetNum) == 0) {
            stop("ERROR: The number of annotated IDs for all functional categories are not from ", minNum, " to ", maxNum, " for the GSEA enrichment method.")
        }

        # collapse rank list
        a <- tapply(geneRankList$score, geneRankList$gene, collapseMethod, na.rm = TRUE)
        geneRankList <- data.frame(gene = names(a), score = as.numeric(a), stringsAsFactors = FALSE)

        gseaRnk <- file.path(projectFolder, paste("Project_", projectName, "_GSEA.rnk", sep = ""))
        write_tsv(geneRankList, gseaRnk, col_names = FALSE)

        outputF <- file.path(projectFolder, paste0("Project_", projectName, "_GSEA/"))
        relativeF <- file.path("./", paste0("Project_", projectName, "_GSEA"))
        if (!dir.exists(outputF) && isOutput) {
            dir.create(outputF)
        }

        inputDf <- prepareInputMatrixGsea(geneRankList, effectiveGeneSet)
        inputDf_list[[i]] <- inputDf
    }
    gseaRes_list <- multiswGsea(inputDf_list,
        thresh_type = "val", perms = perNum,
        min_set_size = minNum, max_set_size = maxNum, p = p,
        nThreads = nThreads, rng_seed = as.integer(format(Sys.time(), "%H%M%S"))
    )
    insig_list <- list()
    sig_list <- list()
    for (j in seq_along(gseaRes_list)) {
        print(paste0("Processing ", j, " ..."))
        gseaRes <- gseaRes_list[[j]]

        # if (saveRawGseaResult) {
        #     saveRDS(gseaRes, file = file.path(outputF, "rawGseaResult.rds"))
        # }
        enrichRes <- gseaRes$Enrichment_Results %>%
            mutate(geneSet = rownames(gseaRes$Enrichment_Results)) %>%
            select(.data$geneSet,
                enrichmentScore = .data$ES, normalizedEnrichmentScore = .data$NES, pValue = .data$p_val, FDR = .data$fdr,
                leadingEdgeNum = .data$leading_edge
            )
        if (sigMethod == "fdr") {
            sig <- filter(enrichRes, .data$FDR < fdrThr)
            insig <- filter(enrichRes, .data$FDR >= fdrThr)
        } else if (sigMethod == "top") {
            enrichRes <- arrange(enrichRes, .data$FDR, .data$pValue)
            if (j == 1) {
                tmpRes <- getTopMetaGseaResults(enrichRes, topThr)
            } else {
                tmpRes <- getTopGseaResults(enrichRes, topThr)
            }
            sig <- tmpRes[[1]]
            insig <- tmpRes[[2]]
        } else {
            warning("WARNING: Invalid significance method ", sigMethod, "!\nDefaulting to top.\n")
            enrichRes <- arrange(enrichRes, .data$FDR, .data$pValue)
            if (j == 1) {
                tmpRes <- getTopMetaGseaResults(enrichRes, topThr)
            } else {
                tmpRes <- getTopGseaResults(enrichRes, topThr)
            }
            sig <- tmpRes[[1]]
            insig <- tmpRes[[2]]
        }
        numSig <- nrow(sig)
        if (numSig == 0) {
            if (sigMethod == "fdr") {
                warning("ERROR: No significant set is identified based on FDR ", fdrThr, "!\n")
            } else {
                warning("ERROR: No significant set is identified based on top ", topThr, "!\n")
            }
            sig_list[[j]] <- NULL
            insig_list[[j]] <- NULL
            next
        }

        if (!is.null(insig)) {
            # insig$leadingEdgeNum <- unname(sapply(insig$geneSet, function(geneSet) {
            # 	rsum <- gseaRes$Running_Sums[, geneSet] # Running sum is a matrix of gene by gene set
            # 	maxPeak <- max(rsum)
            # 	minPeak <- min(rsum)
            # 	if (abs(maxPeak) >= abs(minPeak)) {
            # 		peakIndex <- match(max(rsum), rsum)
            # 		leadingEdgeNum <- sum(gseaRes$Items_in_Set[[geneSet]]$rank <= peakIndex)
            # 	} else {
            # 		peakIndex <- match(min(rsum), rsum)
            # 		leadingEdgeNum <- sum(gseaRes$Items_in_Set[[geneSet]]$rank >= peakIndex)
            # 	}
            # 	return(leadingEdgeNum)
            # }))
        }
        if (j != 1) {
            geneSetName <- geneSetName_list[[j - 1]]
        } else {
            all_gene_set_name <- geneSetName_list[[1]]
            for (i in 2:length(geneSetName_list)) {
                all_gene_set_name <- rbind(all_gene_set_name, geneSetName_list[[i]])
            }
            geneSetName <- all_gene_set_name %>% distinct(.keep_all = TRUE)
        }

        plotSuffix <- ifelse("png" %in% plotFormat, "png", "svg")
        sig <- sig %>%
            left_join(geneSetName, by = "geneSet") %>%
            mutate(size = unname(sapply(geneSet, function(x) nrow(gseaRes$Items_in_Set[[x]]))))
        plot_paths <- list()
        leadingGeneNum <- vector("integer", numSig)
        leadingGenes <- vector("character", numSig)
        if (j == 1) {
            for (i in seq_along(sig$geneSet)) {
                geneSet <- sig[[i, "geneSet"]]
                es <- sig[[i, "enrichmentScore"]]
                genes <- gseaRes$Items_in_Set[[geneSet]] # row name is gene and one column called rank
                leadingGeneNum[[i]] <- nrow(genes)
                leadingGenes[[i]] <- paste(rownames(genes), collapse = ";")
            }
            sig$leadingEdgeNum <- leadingGeneNum
            sig$leadingEdgeId <- leadingGenes
            sig$plotPath <- rep("", length(sig$geneSet))
        } else {
            if (is.null(geneSetDes_list) || length(geneSetDes_list) == 0 || length(geneSetDes_list) < i) {
                geneSetDes <- NULL
            } else {
                geneSetDes <- geneSetDes_list[[i]]
            }
            projectName <- paste0(old_project_name, "_", listNames[j - 1])
            outputF <- file.path(outputDirectory, paste("Project_", old_project_name, "/", projectName, "/plots", sep = ""))
            if (!dir.exists(outputF)) {
                dir.create(outputF)
            }
            for (i in 1:numSig) {
                geneSet <- sig[[i, "geneSet"]]
                es <- sig[[i, "enrichmentScore"]]
                genes <- gseaRes$Items_in_Set[[geneSet]] # row name is gene and one column called rank
                rsum <- gseaRes$Running_Sums[, geneSet]
                peakIndex <- match(ifelse(es > 0, max(rsum), min(rsum)), rsum)
                if (es > 0) {
                    indexes <- genes$rank <= peakIndex
                } else {
                    indexes <- genes$rank >= peakIndex
                }
                leadingGeneNum[[i]] <- sum(indexes)
                leadingGenes[[i]] <- paste(rownames(genes)[indexes], collapse = ";")

                if (isOutput) {
                    # Plot GSEA-like enrichment plot
                    if (!is.null(geneSetDes)) {
                        # same name of variable and column name, use quasi-quotation !!
                        title <- as.character((geneSetDes %>% filter(.data$geneSet == !!geneSet))[1, "description"])
                    } else {
                        title <- geneSet
                    }

                    if (!is.vector(plotFormat)) {
                        plot_paths[[i]] <- plotEnrichmentPlot(title, outputF, geneSet, format = plotFormat, gseaRes$Running_Sums[, geneSet], genes$rank, sortedScores, peakIndex)
                    } else {
                        for (format in plotFormat) {
                            plot_paths[[i]] <- plotEnrichmentPlot(title, outputF, geneSet, format = format, gseaRes$Running_Sums[, geneSet], genes$rank, sortedScores, peakIndex)
                        }
                    }
                }
            }
            sig$plotPath <- unlist(plot_paths)
            sig$leadingEdgeNum <- leadingGeneNum
            sig$leadingEdgeId <- leadingGenes
        }
        sig_list[[j]] <- sig
        insig_list[[j]] <- insig
    }

    return(list(enriched = sig_list, background = insig_list))
}

#' @importFrom svglite svglite
plotEnrichmentPlot <- function(title, outputDir, fileName, format = "png", runningSums, ranks, scores, peakIndex) {
    if (format == "png") {
        output_file <- file.path(outputDir, paste0(sanitizeFileName(fileName), ".png"))
        png(output_file, bg = "transparent", width = 2000, height = 2000)
        cex <- list(main = 5, axis = 2.5, lab = 3.2)
    } else if (format == "svg") {
        output_file <- file.path(outputDir, paste0(sanitizeFileName(fileName), ".svg"))
        svglite(output_file, bg = "transparent", width = 7, height = 7)
        cex <- list(main = 1.5, axis = 0.6, lab = 0.8)
        # svg seems to have a problem with long title (figure margins too large)
        if (!is.na(nchar(title))) {
            if (nchar(title) > 80) {
                title <- paste0(substr(title, 1, 80), "...")
            }
        }
    }
    wrappedTitle <- strwrap(paste0("Enrichment plot: ", title), 60)
    plot.new()
    par(fig = c(0, 1, 0.5, 1), mar = c(0, 6, 6 * length(wrappedTitle), 2), cex.axis = cex$axis, cex.main = cex$main, cex.lab = cex$lab, lwd = 2, new = TRUE)
    plot(1:length(runningSums), runningSums,
        type = "l", main = paste(wrappedTitle, collapse = "\n"),
        xlab = "", ylab = "Enrichment Score", xaxt = "n", lwd = 3
    )
    abline(v = peakIndex, lty = 3)
    par(fig = c(0, 1, 0.35, 0.5), mar = c(0, 6, 0, 2), new = TRUE)
    plot(ranks, rep(1, length(ranks)),
        type = "h",
        xlim = c(1, length(scores)), ylim = c(0, 1), axes = FALSE, ann = FALSE
    )
    par(fig = c(0, 1, 0, 0.35), mar = c(6, 6, 0, 2), cex.axis = cex$axis, cex.lab = cex$lab, new = TRUE)
    # use polygon to greatly reduce file size of SVG
    plot(1:length(scores), scores,
        type = "n",
        ylab = "Ranked list metric", xlab = "Rank in Ordered Dataset"
    )
    polygon(c(1, 1:length(scores), length(scores)), c(0, scores, 0), col = "black")
    abline(v = peakIndex, lty = 3)
    dev.off()
    return(output_file)
}
