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
            allEnrichedSig <- enrichedSig
            repAdded <- FALSE
            if (organism != "others") {
                if (!is.null(enrichedSig) && reportNum < nrow(enrichedSig)) {
                    if (enrichMethod == "ORA") {
                        enrichedSig <- enrichedSig[1:reportNum, ]
                    } else if (enrichMethod == "GSEA") {
                        enrichedSig <- getTopGseaResults(enrichedSig, reportNum / 2)[[1]]
                    }
                    # Add representatives if they are not in top ReportNum. So could be more if ReportNum.is small and high redundancy in top
                    numRes <- nrow(enrichedSig)
                    enrichedSig <- keepRep(enrichedSig, allEnrichedSig, clusters$ap$representatives)
                    enrichedSig <- keepRep(enrichedSig, allEnrichedSig, clusters$wsc$representatives)
                    repAdded <- nrow(enrichedSig) > numRes
                }
                standardId <- interestingGeneMap$standardId
                if (enrichMethod == "ORA") {
                    interestGeneList <- unique(interestingGeneMap$mapped[[standardId]])
                    numAnnoRefUserId <- length(intersect(
                        interestGeneList,
                        intersect(referenceGeneList, geneSet$gene)
                    ))
                }

                ##### Summary Tab ########
                bodyContent <- summaryDescription(projectName, organism, interestGeneFile, interestGene, interestGeneType, enrichMethod, enrichDatabase, enrichDatabaseFile, enrichDatabaseType, enrichDatabaseDescriptionFile, interestingGeneMap, referenceGeneList, referenceGeneFile, referenceGene, referenceGeneType, referenceSet, minNum, maxNum, sigMethod, fdrThr, topThr, fdrMethod, allEnrichedSig, reportNum, perNum, p, geneSet, repAdded, numAnnoRefUserId, hostName)

                ########### GOSlim summary #########################
                if (standardId == "entrezgene") {
                    bodyContent <- paste(bodyContent, goSlimReport(projectName), sep = "\n")
                }

                ############ Enrichment result ##################
                if (!is.null(enrichedSig)) {
                    bodyContent <- paste(bodyContent, enrichResultSection(enrichMethod, enrichedSig, geneSet, geneSetDes, geneSetDag, geneSetNet, clusters), seq = "\n")
                    if (!is.null(geneSetDag)) {
                        if (!is.vector(geneSetDag)) {
                            # for backward compatibility, it is unlisted for single dataset
                            geneSetDag <- list(geneSetDag)
                            names(geneSetDag) <- ifelse(is.character(enrichDatabase), enrichDatabase, gsub(".gmt", "", basename(enrichDatabaseFile), fixed = TRUE))
                        }
                        for (name in names(geneSetDag)) {
                            dag <- geneSetDag[[name]]
                            if (is.null(dag)) {
                                # dagJson[[name]] <- list(NULL)
                                next
                            }
                            dagRes <- expandDag(enrichedSig$geneSet, dag)
                            dagEdges <- dagRes$edges
                            dagNodes <- getDagNodes(enrichedSig, dagRes$allNodes, geneSetDes, enrichMethod, dagColor)
                            dagJson[[name]] <- c(dagEdges, dagNodes)
                        }
                    }
                }
            } else {
                ########### Organism is others. No mapping information #############
                ############# Summary for the analysis ###################
                if (enrichMethod == "ORA") {
                    numAnnoRefUserId <- length(intersect(
                        interestingGeneMap,
                        intersect(referenceGeneList, geneSet$gene)
                    ))
                }
                if (!is.null(enrichedSig) && reportNum < nrow(enrichedSig)) {
                    if (enrichMethod == "ORA") {
                        enrichedSig <- enrichedSig[1:reportNum, ]
                    } else if (enrichMethod == "GSEA") {
                        enrichedSig <- getTopGseaResults(enrichedSig, reportNum / 2)[[1]]
                    }
                    # Add representatives if they are not in top ReportNum. So could be more if ReportNum.is small and high redundancy in top
                    numRes <- nrow(enrichedSig)
                    enrichedSig <- keepRep(enrichedSig, allEnrichedSig, clusters$ap$representatives)
                    enrichedSig <- keepRep(enrichedSig, allEnrichedSig, clusters$wsc$representatives)
                    repAdded <- nrow(enrichedSig) > numRes
                }

                bodyContent <- summaryDescription(projectName, organism, interestGeneFile, interestGene, interestGeneType, enrichMethod, enrichDatabase, enrichDatabaseFile, enrichDatabaseType, enrichDatabaseDescriptionFile, interestingGeneMap, referenceGeneList, referenceGeneFile, referenceGene, referenceGeneType, referenceSet, minNum, maxNum, sigMethod, fdrThr, topThr, fdrMethod, allEnrichedSig, reportNum, perNum, p, geneSet, repAdded, numAnnoRefUserId, hostName)

                ############## Enrich Result ################
                if (!is.null(enrichedSig)) {
                    bodyContent <- paste(bodyContent, enrichResultSection(enrichMethod, enrichedSig, geneSet, geneSetDes, geneSetDag, geneSetNet, clusters), seq = "\n")
                }
                standardId <- NULL
            }
            if (is.null(enrichedSig)) {
                enrichedSig <- data.frame()
            }
            if (is.null(background)) {
                background <- data.frame()
            }
            version <- packageVersion("WebGestaltR")
            # use major and minor version numbers for JS lib. If API changes, version should be bumped
            # patch number should not matter
            version <- paste(version[1, 1], version[1, 2], sep = ".")
            hasGeneSetDag <- !is.null(geneSetDag)
            hasCytoscape <- hasGeneSetDag || !is.null(geneSetNet) # DAG or network needs cytoscape
            allDbNames <- unlist(c(enrichDatabase, unname(sapply(enrichDatabaseFile, function(x) {
                gsub(".gmt", "", basename(x), fixed = TRUE)
            })))) # sapply on NULL will return a list
            tabs[[i]] <- list(title = listNames[i], bodyContent = bodyContent)
        }
    }

    allContent <- "<b-tabs>\n"
    for (i in seq_along(tabs)) {
      allContent <- paste(allContent, "<b-tab title=\"", tabs[[i]]$title, "\" active>", tabs[[i]]$bodyContent, "</b-tab>\n", sep = "")
    }
    allContent <- paste(allContent, "</b-tabs>\n", sep = "")
    # tabs <- list(tabs = tabs)
    print(toJSON(tabs))
    header <- readLines(system.file("templates/header.mustache", package = "WebGestaltR"))
    footer <- readLines(system.file("templates/footer.mustache", package = "WebGestaltR"))
    template <- readLines(system.file("templates/meta_template.mustache", package = "WebGestaltR"))
    data <- list(
        hostName = hostName, allContent = allContent,
        organism = organism, enrichDatabaseJson = toJSON(allDbNames, auto_unbox = TRUE),
        sigJson = toJSON(enrichedSig, digits = 16), insigJson = toJSON(background, digits = 16),
        dagJson = toJSON(dagJson, auto_unbox = TRUE), hasGeneSetDag = hasGeneSetDag, version = version,
        clusterJson = toJSON(clusters), hasCytoscape = hasCytoscape,
        geneTableJson = toJSON(geneTables), standardId = standardId, numAnnoRefUserId = numAnnoRefUserId,
        methodIsGsea = enrichMethod == "GSEA", hasGeneSetDes = !is.null(geneSetDes)
    )
    cat(whisker.render(template, data, partials = list(header = header, footer = footer)), file = outputHtmlFile)
}
