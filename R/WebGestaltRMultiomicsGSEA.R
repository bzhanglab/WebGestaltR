#' @title Multi-omics GSEA
#' importFrom dplyr bind_rows left_join arrange select desc %>% mutate distinct
#' importFrom readr write_tsv
#' @inheritParams WebGestaltRMultiOmics
WebGestaltRMultiOmicsGSEA <- function(analyteLists = NULL, analyteListFiles = NULL, analyteTypes = NULL, enrichMethod = "GSEA", organism = "hsapiens",
                                      enrichDatabase = NULL, enrichDatabaseFile = NULL, enrichDatabaseType = NULL, enrichDatabaseDescriptionFile = NULL,
                                      collapseMethod = "mean", minNum = 10, maxNum = 500, fdrMethod = "BH", sigMethod = "fdr", fdrThr = 0.05,
                                      topThr = 10, reportNum = 100, setCoverNum = 10, perNum = 1000, gseaP = 1, isOutput = TRUE, outputDirectory = getwd(),
                                      projectName = NULL, dagColor = "binary", saveRawGseaResult = FALSE, gseaPlotFormat = "png", nThreads = 1, cache = NULL,
                                      hostName = "https://www.webgestalt.org/", useWeightedSetCover = TRUE, useAffinityPropagation = FALSE,
                                      usekMedoid = FALSE, kMedoid_k = 25, isMetaAnalysis = TRUE, mergeMethod = "mean", normalizationMethod = "rank",
                                      listNames = NULL) {
    projectDir <- file.path(outputDirectory, paste0("Project_", projectName))
    cat("Performing multi-omics GSEA\nLoading the functional categories...\n")
    all_sets <- .load_meta_gmt(enrichDatabase, enrichDatabaseFile, enrichDatabaseDescriptionFile, enrichDatabaseType, analyteLists, analyteListFiles, analyteTypes, organism, cache, hostName)
    cat("Loading the ID lists...\n")
    interest_lists <- list()
    interestGeneMaps <- list()
    if (is.null(analyteLists)) {
        for (i in seq_along(analyteListFiles)) {
            interestingGeneMap <- loadInterestGene(
                organism = organism, dataType = "rnk", inputGeneFile = analyteListFiles[[i]], inputGene = NULL, geneType = analyteTypes[[i]],
                collapseMethod = collapseMethod, cache = cache, hostName = hostName, geneSet = all_sets[["geneSet"]][[i]]
            )
            interestGeneMaps[[i]] <- interestingGeneMap
            if (organism == "others") {
                interestGeneList <- unique(interestingGeneMap)
                interest_lists[[i]] <- interestGeneList
            } else {
                interestStandardId <- interestingGeneMap$standardId
                interestGeneList <- unique(interestingGeneMap$mapped[[interestStandardId]])
                interest_lists[[i]] <- interestGeneList
            }
        }
    } else {
        for (i in seq_along(analyteLists)) {
            interestingGeneMap <- loadInterestGene(
                organism = organism, dataType = "rnk", inputGeneFile = NULL, inputGene = analyteLists[[i]], geneType = analyteTypes[[i]],
                collapseMethod = collapseMethod, cache = cache, hostName = hostName, geneSet = all_sets[["geneSet"]][[i]]
            )
            interestGeneMaps[[i]] <- interestingGeneMap
            if (organism == "others") {
                interestGeneList <- unique(interestingGeneMap)
                interest_lists[[i]] <- interestGeneList
            } else {
                interestStandardId <- interestingGeneMap$standardId
                interestGeneList <- unique(interestingGeneMap$mapped[[interestStandardId]])
                interest_lists[[i]] <- interestGeneList
            }
        }
    }

    cat("Running multi-omics GSEA...\n")

    multiGseaEnrichment(hostName, outputDirectory, projectName, interest_lists, all_sets[["geneSet"]], all_sets[["geneSetDes"]],
        collapseMethod = "mean",
        minNum = 10, maxNum = 500, sigMethod = "fdr", fdrThr = 0.05, topThr = 10, perNum = 1000, p = 1, isOutput = TRUE,
        saveRawGseaResult = FALSE, plotFormat = "png", nThreads = 1, listNames = NULL
    )

    print("Ran successfully!")
}
