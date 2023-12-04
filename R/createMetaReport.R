#' createMetaReport
#'
#' Generate HTML report for ORA and GSEA MetaAnalysis
#'
#' @importFrom jsonlite toJSON
#' @importFrom whisker whisker.render
#'
#' @keywords internal
#'
createMetaReport <- function(hostName, outputDirectory, organism = "hsapiens", projectName, enrichMethod, geneSet_list, geneSetDes_list,
                             geneSetDag_list, geneSetNet_list, interestingGeneMap_list, referenceGeneList_list, enrichedSig_list, geneTables_list, clusters_list,
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
    for (i in seq_along(interestingGeneMap_list)) {
        if (i == 0) {

        } else {
            print(paste("Processing dataset", i, "of", length(interestingGeneMap_list)))
            enrichedSig <- enrichedSig_list[[i + 1]]
            geneSet <- geneSet_list[[i]]
            if (!is.null(geneSetDes_list) && length(geneSetDes_list) > 0) {
                geneSetDes <- geneSetDes_list[[i]]
            } else {
                geneSetDes <- NULL
            }
            if (!is.null(geneSetDag_list) && length(geneSetDag_list) > 0) {
                geneSetDag <- geneSetDag_list[[i]]
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
            geneTables <- geneTables_list[[i]]
            clusters <- clusters_list[[i]]
            background <- background_list[[i]]
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
            numAnnoRefUserId <- NULL
            dagJson <- list()
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
                topThr = topThr, reportNum = reportNum, dagColor = dagColor, outputHtmlFile = partial_output, is_meta = TRUE
            )
            tabs[[i]] <- list(title = listNames[[i]], path = relative_path)
        }
    }

    allContent <- "<b-tabs>\n"
    for (i in seq_along(tabs)) {
        allContent <- paste(allContent, "<b-tab-item label=\"", tabs[[i]]$title, "\"><iframe style=\"width: 100%\" src=\"", tabs[[i]]$path, "\"></iframe></b-tab-item>\n", sep = "")
    }
    allContent <- paste(allContent, "</b-tabs>\n", sep = "")
    # tabs <- list(tabs = tabs)
    print(toJSON(tabs))
    header <- readLines(system.file("templates/header.mustache", package = "WebGestaltR"))
    footer <- readLines(system.file("templates/footer.mustache", package = "WebGestaltR"))
    template <- readLines(system.file("templates/meta_template.mustache", package = "WebGestaltR"))
    data <- list(
        hostName = hostName, allContent = allContent, version = version
    )
    cat(whisker.render(template, data, partials = list(header = header, footer = footer)), file = outputHtmlFile)
}
