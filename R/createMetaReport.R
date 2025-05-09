#' createMetaReport
#'
#' Generate HTML report for ORA and GSEA MetaAnalysis
#'
#' @importFrom jsonlite toJSON
#' @importFrom whisker whisker.render
#'
#' @keywords internal
#'
createMetaReport <- function(hostName = NULL, outputDirectory = NULL, organism = "hsapiens", projectName = NULL, enrichMethod = NULL, geneSet_list = NULL, geneSetDes_list = NULL,
                             geneSetDag_list = NULL, geneSetNet_list = NULL, interestingGeneMap_list = NULL, referenceGeneList_list = NULL, enrichedSig_list = NULL, geneTables_list = NULL, clusters_list = NULL,
                             background_list, enrichDatabase_list = NULL, enrichDatabaseFile_list = NULL, enrichDatabaseType_list = NULL,
                             enrichDatabaseDescriptionFile_list = NULL, interestGeneFile_list = NULL, interestGene_list = NULL,
                             interestGeneType_list = NULL, collapseMethod = "mean", referenceGeneFile_list = NULL, referenceGene_list = NULL,
                             referenceGeneType_list = NULL, referenceSet_list = NULL, minNum = 10, maxNum = 500, fdrMethod = "BH", sigMethod = "fdr", fdrThr = 0.05,
                             topThr = 10, reportNum = 20, perNum = 1000, p = 1, dagColor = "binary", listNames = NULL) {
    outputHtmlFile <- file.path(outputDirectory, paste0("Project_", projectName), paste0("Report_", projectName, ".html"))
    version <- packageVersion("WebGestaltR")
    # use major and minor version numbers for JS lib. If API changes, version should be bumped
    # patch number should not matter
    version <- paste(version[1, 1], version[1, 2], sep = ".")
    # if hostname starts with "file://", it is used as WebGestaltReporter
    if (startsWith(hostName, "file://")) {
        # change back hostName for web assets and browsers will cache it.
        hostName <- "https://www.webgestalt.org"
    }
    tabs <- c()
    for (j in seq_along(enrichedSig_list)) {
        if (j == 1) {
            print(paste("Generating report for dataset", j, "of", length(enrichedSig_list)))
            enrichedSig <- enrichedSig_list[[j]]
            geneSet <- NULL
            if (!is.null(geneSetDes_list) && length(geneSetDes_list) > 0) {
                geneSetDes <- NULL
            } else {
                geneSetDes <- NULL
            }
            if (!is.null(geneSetDag_list) && length(geneSetDag_list) > 0) {
                geneSetDag <- NULL
            } else {
                geneSetDag <- NULL
            }
            if (!is.null(geneSetNet_list) && length(geneSetNet_list) > 0) {
                geneSetNet <- NULL
            } else {
                geneSetNet <- NULL
            }
            interestingGeneMap <- NULL
            referenceGeneList <- NULL
            geneTables <- geneTables_list[[j]]
            clusters <- clusters_list[[j]]
            background <- background_list[[j]]
            if (length(enrichDatabase_list) == 1) {
                enrichDatabase <- enrichDatabase_list
            } else {
                enrichDatabase <- NULL
            }
            if (length(enrichDatabaseFile_list) == 1) {
                enrichDatabaseFile <- enrichDatabaseFile_list
            } else {
                enrichDatabaseFile <- NULL
            }
            if (length(enrichDatabaseType_list) == 1) {
                enrichDatabaseType <- enrichDatabaseType_list
            } else {
                enrichDatabaseType <- NULL
            }
            if (length(enrichDatabaseDescriptionFile_list) == 1) {
                enrichDatabaseDescriptionFile <- enrichDatabaseDescriptionFile_list
            } else {
                enrichDatabaseDescriptionFile <- NULL
            }
            interestGeneFile <- NULL
            interestGene <- NULL
            interestGeneType <- NULL
            referenceGeneFile <- NULL
            referenceGene <- NULL
            referenceGeneType <- NULL
            referenceSet <- NULL
            relative_path <- paste0("./Report_", projectName, "_meta", ".html")
            partial_output <- file.path(outputDirectory, paste0("Project_", projectName), paste0("Report_", projectName, "_meta", ".html"))
            pvals <- unlist(enrichedSig$pValue)
            nes <- unlist(enrichedSig$normalizedEnrichmentScore)
            logp <- c()
            safe_sign <- function(x) {
                if (x > 0) {
                    return(1)
                } else if (x < 0) {
                    return(-1)
                } else {
                    return(1)
                }
            }
            for (i in 1:length(pvals)) {
                x <- pvals[i]
                if (enrichMethod == "ORA") {
                    if (abs(x) <= 10 * .Machine$double.eps) {
                        logp[i] <- 16
                    } else {
                        logp[i] <- -log10(abs(x))
                    }
                } else if (enrichMethod == "GSEA") {
                    if (sign(nes[i]) == 0) {
                        nes[i] <- 1
                    }
                    if (abs(x) <= 2 * .Machine$double.eps) {
                        logp[i] <- -log10(.Machine$double.eps) * safe_sign(nes[i])
                    } else {
                        logp[i] <- -log10(abs(x)) * safe_sign(nes[i])
                    }
                }
                if (is.na(logp[i])) {
                    logp[i] <- .Machine$double.eps
                } else if (is.infinite(logp[i])) {
                    logp[i] <- .Machine$double.eps
                }
            }
            enrichedSig$logp <- unlist(logp)

            createMetaSummaryReport(
                hostName = hostName, outputDirectory = outputDirectory, organism = organism,
                projectName = projectName, enrichMethod = enrichMethod, geneSet = geneSet,
                geneSetDes = geneSetDes, geneSetDag = geneSetDag, geneSetNet = geneSetNet,
                interestingGeneMap = interestingGeneMap, referenceGeneList = referenceGeneList,
                enrichedSig = enrichedSig, background = background, geneTables = geneTables,
                clusters = clusters, enrichDatabase = enrichDatabase,
                enrichDatabaseFile = enrichDatabaseFile, enrichDatabaseType = enrichDatabaseType,
                enrichDatabaseDescriptionFile = enrichDatabaseDescriptionFile,
                interestGeneFile = interestGeneFile, interestGene = interestGene,
                interestGeneType = interestGeneType, collapseMethod = collapseMethod,
                referenceGeneFile = referenceGeneFile, referenceGene = referenceGene,
                referenceGeneType = referenceGeneType, referenceSet = referenceSet, minNum = minNum,
                maxNum = maxNum, fdrMethod = fdrMethod, sigMethod = sigMethod, fdrThr = fdrThr,
                topThr = topThr, reportNum = reportNum, dagColor = dagColor, outputHtmlFile = partial_output, is_meta = TRUE, listNames = listNames
            )
            tabs[[j]] <- list(title = "All", path = relative_path)
        } else {
            i <- j - 1 # offset for lists that don't have meta-analysis
            print(paste("Generating report for dataset", j, "of", length(enrichedSig_list)))
            enrichedSig <- enrichedSig_list[[j]]
            geneSet <- geneSet_list[[i]]
            if (!is.null(geneSetDes_list) && length(geneSetDes_list) > 0 && i <= length(geneSetDes_list)) {
                geneSetDes <- geneSetDes_list[[i]]
            } else {
                geneSetDes <- NULL
            }
            if (!is.null(geneSetDag_list) && length(geneSetDag_list) > 0) {
                #  geneSetDag <- geneSetDag_list[[i]]
                geneSetDag <- NULL
            } else {
                geneSetDag <- NULL
            }
            if (!is.null(geneSetNet_list) && length(geneSetNet_list) > 0) {
                geneSetNet <- geneSetNet_list[[i]]
            } else {
                geneSetNet <- NULL
            }
            interestingGeneMap <- interestingGeneMap_list[[i]]
            referenceGeneList <- referenceGeneList_list[[i]]
            geneTables <- geneTables_list[[j]]
            clusters <- clusters_list[[j]]
            background <- background_list[[j]]
            if (length(enrichDatabase_list) == 1) {
                enrichDatabase <- enrichDatabase_list
            } else {
                enrichDatabase <- enrichDatabase_list[[i]]
            }
            if (length(enrichDatabaseFile_list) == 1) {
                enrichDatabaseFile <- enrichDatabaseFile_list
            } else {
                enrichDatabaseFile <- enrichDatabaseFile_list[[i]]
            }
            if (length(enrichDatabaseType_list) == 1) {
                enrichDatabaseType <- enrichDatabaseType_list
            } else {
                enrichDatabaseType <- enrichDatabaseType_list[[i]]
            }
            if (length(enrichDatabaseDescriptionFile_list) == 1) {
                enrichDatabaseDescriptionFile <- enrichDatabaseDescriptionFile_list
            } else {
                enrichDatabaseDescriptionFile <- enrichDatabaseDescriptionFile_list[[i]]
            }
            interestGeneFile <- interestGeneFile_list[[i]]
            interestGene <- interestGene_list[[i]]
            interestGeneType <- interestGeneType_list[[i]]
            referenceGeneFile <- referenceGeneFile_list[[i]]
            referenceGene <- referenceGene_list[[i]]
            referenceGeneType <- referenceGeneType_list[[i]]
            referenceSet <- referenceSet_list[[i]]
            relative_path <- paste0("./Report_", projectName, "_", i, ".html")
            partial_output <- file.path(outputDirectory, paste0("Project_", projectName), paste0("Report_", projectName, "_", i, ".html"))
            createReport(
                hostName = hostName, outputDirectory = outputDirectory, organism = organism,
                projectName = projectName, enrichMethod = enrichMethod, geneSet = geneSet,
                geneSetDes = geneSetDes, geneSetDag = geneSetDag, geneSetNet = geneSetNet,
                interestingGeneMap = interestingGeneMap, referenceGeneList = referenceGeneList,
                enrichedSig = enrichedSig, background = background, geneTables = geneTables,
                clusters = clusters, enrichDatabase = enrichDatabase,
                enrichDatabaseFile = enrichDatabaseFile, enrichDatabaseType = enrichDatabaseType,
                enrichDatabaseDescriptionFile = enrichDatabaseDescriptionFile,
                interestGeneFile = interestGeneFile, interestGene = interestGene,
                interestGeneType = interestGeneType, collapseMethod = collapseMethod,
                referenceGeneFile = referenceGeneFile, referenceGene = referenceGene,
                referenceGeneType = referenceGeneType, referenceSet = referenceSet, minNum = minNum,
                maxNum = maxNum, fdrMethod = fdrMethod, sigMethod = sigMethod, fdrThr = fdrThr,
                topThr = topThr, reportNum = reportNum, dagColor = dagColor, outputHtmlFile = partial_output, listName= listNames[[i]],is_meta = TRUE
            )
            tabs[[j]] <- list(title = listNames[[i]], path = relative_path)
        }
    }

    allContent <- "<b-tabs>\n"
    for (i in seq_along(tabs)) {
        allContent <- paste(allContent, "<b-tab-item label=\"", tabs[[i]]$title, "\"><iframe style=\"width: 100%\" src=\"", tabs[[i]]$path, "\"></iframe></b-tab-item>\n", sep = "")
    }
    allContent <- paste(allContent, "</b-tabs>\n", sep = "")
    header <- readLines(system.file("templates/header.mustache", package = "WebGestaltR"))
    footer <- readLines(system.file("templates/footer.mustache", package = "WebGestaltR"))
    template <- readLines(system.file("templates/meta_template.mustache", package = "WebGestaltR"))
    data <- list(
        hostName = hostName, allContent = allContent, version = version, html_title = paste(listNames, collapse = "+")
    )
    cat(whisker.render(template, data, partials = list(header = header, footer = footer)), file = outputHtmlFile)
}
