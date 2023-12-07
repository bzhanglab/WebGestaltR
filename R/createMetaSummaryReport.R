#' createMetaSummaryReport
#'
#' Generate HTML report for ORA and GSEA
#'
#' @importFrom jsonlite toJSON
#' @importFrom whisker whisker.render
#'
#' @keywords internal
#'
createMetaSummaryReport <- function(hostName, outputDirectory, organism = "hsapiens", projectName, enrichMethod, geneSet, geneSetDes, geneSetDag, geneSetNet,
                                    interestingGeneMap, referenceGeneList, enrichedSig, geneTables, clusters, background, enrichDatabase = NULL, enrichDatabaseFile = NULL, enrichDatabaseType = NULL,
                                    enrichDatabaseDescriptionFile = NULL, interestGeneFile = NULL, interestGene = NULL, interestGeneType = NULL, collapseMethod = "mean", referenceGeneFile = NULL,
                                    referenceGene = NULL, referenceGeneType = NULL, referenceSet = NULL, minNum = 10, maxNum = 500, fdrMethod = "BH", sigMethod = "fdr", fdrThr = 0.05, topThr = 10,
                                    reportNum = 20, perNum = 1000, p = 1, dagColor = "binary", is_meta = TRUE, outputHtmlFile = NULL, listNames = NULL) {
    if (is.null(outputHtmlFile)) {
        outputHtmlFile <- file.path(outputDirectory, paste0("Project_", projectName), paste0("Report_", projectName, ".html"))
    }
    # if hostname starts with "file://", it is used as WebGestaltReporter
    if (startsWith(hostName, "file://")) {
        # change back hostName for web assets and browsers will cache it.
        hostName <- "https://www.webgestalt.org"
    }

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
        bodyContent <- metaSummaryDescription(
            projectName, organism, interestGeneFile, interestGene, interestGeneType, enrichMethod, enrichDatabase,
            enrichDatabaseFile, enrichDatabaseType, enrichDatabaseDescriptionFile, interestingGeneMap, referenceGeneList,
            referenceGeneFile, referenceGene, referenceGeneType, referenceSet, minNum, maxNum, sigMethod, fdrThr, topThr,
            fdrMethod, allEnrichedSig, reportNum, perNum, p, geneSet, repAdded, numAnnoRefUserId, hostName, listNames
        )

        ############ Enrichment result ##################
        if (!is.null(enrichedSig)) {
            bodyContent <- paste(bodyContent, metaEnrichResultSection(enrichMethod, enrichedSig, geneSet, geneSetDes, geneSetDag, geneSetNet, clusters), seq = "\n")
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

        bodyContent <- metaSummaryDescription(
            projectName, organism, interestGeneFile, interestGene, interestGeneType, enrichMethod, enrichDatabase,
            enrichDatabaseFile, enrichDatabaseType, enrichDatabaseDescriptionFile, interestingGeneMap, referenceGeneList,
            referenceGeneFile, referenceGene, referenceGeneType, referenceSet, minNum, maxNum, sigMethod, fdrThr, topThr,
            fdrMethod, allEnrichedSig, reportNum, perNum, p, geneSet, repAdded, numAnnoRefUserId, hostName, listNames
        )

        ############## Enrich Result ################
        if (!is.null(enrichedSig)) {
            bodyContent <- paste(bodyContent, metaEnrichResultSection(enrichMethod, enrichedSig, geneSet, geneSetDes, geneSetDag, geneSetNet, clusters), seq = "\n")
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
    hasGeneSetDag <- FALSE
    hasCytoscape <- hasGeneSetDag || !is.null(geneSetNet) # DAG or network needs cytoscape
    allDbNames <- unlist(c(enrichDatabase, unname(sapply(enrichDatabaseFile, function(x) {
        gsub(".gmt", "", basename(x), fixed = TRUE)
    })))) # sapply on NULL will return a list

    template <- readLines(system.file("templates/meta_analysis_template.mustache", package = "WebGestaltR"))
    data <- list(
        hostName = hostName, bodyContent = bodyContent,
        organism = organism, enrichDatabaseJson = toJSON(allDbNames, auto_unbox = TRUE),
        sigJson = toJSON(enrichedSig, digits = 16), insigJson = toJSON(background, digits = 16),
        dagJson = toJSON(dagJson, auto_unbox = TRUE), hasGeneSetDag = FALSE, version = version,
        clusterJson = toJSON(clusters), hasCytoscape = hasCytoscape,
        geneTableJson = toJSON(geneTables), standardId = standardId, numAnnoRefUserId = numAnnoRefUserId,
        methodIsGsea = enrichMethod == "GSEA", hasGeneSetDes = !is.null(geneSetDes), listDataJson = toJSON(listNames)
    )
    cat(whisker.render(template, data), file = outputHtmlFile)
}
