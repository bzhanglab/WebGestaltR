#' @title Multi-omics GSEA
#' importFrom dplyr bind_rows left_join arrange select desc %>% mutate distinct
#' importFrom readr write_tsv
#' @inheritParams WebGestaltRMultiOmics
WebGestaltRMultiOmicsGSEA <- function(analyteLists = NULL, analyteListFiles = NULL, analyteTypes = NULL, enrichMethod = "GSEA", organism = "hsapiens",
                                      enrichDatabase = NULL, enrichDatabaseFile = NULL, enrichDatabaseType = NULL, enrichDatabaseDescriptionFile = NULL,
                                      collapseMethod = "mean", minNum = 10, maxNum = 500, fdrMethod = "BH", sigMethod = "fdr", fdrThr = 0.05,
                                      topThr = 10, reportNum = 100, setCoverNum = 10, perNum = 1000, gseaP = 1, isOutput = TRUE, outputDirectory = getwd(),
                                      projectName = NULL, dagColor = "binary", saveRawGseaResult = FALSE, gseaPlotFormat = "png", nThreads = 1,
                                      cache = NULL, hostName = "https://www.webgestalt.org/", useWeightedSetCover = TRUE, useAffinityPropagation = FALSE,
                                      usekMedoid = FALSE, kMedoid_k = 25, isMetaAnalysis = TRUE, mergeMethod = "mean", normalizationMethod = "rank",
                                      listNames = NULL) {
    projectDir <- file.path(outputDirectory, paste0("Project_", projectName))
    cat("Performing multi-omics GSEA\nLoading the functional categories...\n")
    all_sets <- .load_meta_gmt(enrichDatabase, enrichDatabaseFile, enrichDatabaseDescriptionFile, enrichDatabaseType, analyteLists, analyteListFiles, analyteTypes,
                               organism, cache, hostName)
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
            } else {
                interestStandardId <- interestingGeneMap$standardId
                interestGeneList <- interestingGeneMap$mapped %>% select(interestStandardId, .data$score) %>% distinct()
            }
            interest_lists[[i]] <- interestGeneList
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
            } else {
                interestStandardId <- interestingGeneMap$standardId
                interestGeneList <- interestingGeneMap$mapped %>% select(interestStandardId, .data$score) %>% distinct()
            }
            interest_lists[[i]] <- interestGeneList
        }
    }

    cat("Running multi-omics GSEA...\n")

    multiGseaEnrichment(hostName = hostName, outputDirectory = outputDirectory, projectName = projectName, geneRankList_list = interest_lists, geneSet_list = all_sets[["geneSet"]], geneSetDes_list = all_sets[["geneSetDes"]],
                        collapseMethod = "mean", minNum = minNum, maxNum = maxNum, sigMethod = sigMethod, fdrThr = fdrThr, topThr = topThr, perNum = perNum, p = p,
                        isOutput = isOutput, saveRawGseaResult = saveRawGseaResult, plotFormat = gseaPlotFormat, nThreads = nThreads, listNames = listNames
    )

    print("Ran successfully!")
}
