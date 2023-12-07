#' metaSummaryDescription
#'
#' Render job summary section
#'
#' @importFrom whisker whisker.render
#' @keywords internal
#'
metaSummaryDescription <- function(projectName, organism, interestGeneFile, interestGene, interestGeneType, enrichMethod, enrichDatabase,
                                   enrichDatabaseFile, enrichDatabaseType, enrichDatabaseDescriptionFile, interestingGeneMap, referenceGeneList,
                                   referenceGeneFile, referenceGene, referenceGeneType, referenceSet, minNum, maxNum, sigMethod, fdrThr, topThr,
                                   fdrMethod, enrichedSig, reportNum, perNum, p, geneSet, repAdded, numAnnoRefUserId, hostName, listNames) {
    template <- readLines(system.file("templates/metaSummary.mustache", package = "WebGestaltR"))
    enrichDatabaseInfo <- list()
    for (enrichDb in enrichDatabase) {
        if (!is.null(enrichDb)) {
            enrichDatabaseInfo <- c(enrichDatabaseInfo, list(list(isBuiltIn = TRUE, enrichDatabase = enrichDb)))
        }
    }
    for (i in seq_along(enrichDatabaseFile)) {
        if (!is.null(enrichDatabaseFile[[i]])) {
            enrichDatabaseInfo <- c(enrichDatabaseInfo, list(list(
                isBuiltIn = FALSE, enrichDatabaseFile = enrichDatabaseFile[[i]],
                hasEnrichDatabaseDescriptionFile = is.null(enrichDatabaseDescriptionFile[[i]]),
                enrichDatabaseDescriptionFile = enrichDatabaseDescriptionFile[[i]]
            )))
        }
    }
    if (organism != "others") {
        standardId <- unname(interestingGeneMap$standardId)
        data <- list(
            projectName = projectName, enrichMethod = enrichMethod, organism = organism, organismIsOthers = FALSE,
            enrichDatabaseInfo = enrichDatabaseInfo, enrichDatabaseType = enrichDatabaseType,
            hasInterestGeneFile = !is.null(interestGeneFile),
            interestGeneFileBase = ifelse(is.null(interestGeneFile), "", basename(interestGeneFile)), interestGeneType = interestGeneType,
            idIsEntrezGene = standardId == "entrezgene", hostName = hostName
        )
    } else {
        data <- list(
            projectName = projectName, enrichMethod = enrichMethod, organism = organism, organismIsOthers = TRUE,
            enrichDatabaseInfo = enrichDatabaseInfo, enrichDatabaseType = enrichDatabaseType,
            hasInterestGeneFile = !is.null(interestGeneFile), interestGeneFileBase = ifelse(is.null(interestGeneFile), "",
                basename(interestGeneFile)
            ), interestGeneType = interestGeneType, idIsEntrezGene = FALSE
        )
    }
    data["listCount"] <- length(listNames)
    return(whisker.render(template, data))
}
