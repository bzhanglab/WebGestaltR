WebGestaltRMultiOmics <- function(analyteLists = NULL, analyteListFiles = NULL, analyteTypes = NULL, enrichMethod = "ORA", organism = "hsapiens", enrichDatabase = NULL,
                                  enrichDatabaseFile = NULL, enrichDatabaseType = NULL, enrichDatabaseDescriptionFile = NULL, interestGeneFile = NULL,
                                  collapseMethod = "mean", minNum = 10, maxNum = 500, fdrMethod = "BH", sigMethod = "fdr", fdrThr = 0.05,
                                  topThr = 10, reportNum = 20, setCoverNum = 10, perNum = 1000, p = 1, isOutput = TRUE, outputDirectory = getwd(),
                                  projectName = NULL, dagColor = "binary", saveRawGseaResult = FALSE, plotFormat = "png", nThreads = 1, cache = NULL,
                                  hostName = "https://www.webgestalt.org/", useWeightedSetCover = TRUE, useAffinityPropagation = FALSE,
                                  usekMedoid = FALSE, kMedoid_k = 10, isMetaAnalysis = TRUE, mergeMethod = "mean", normalizationMethod = "rank",
                                  referenceLists = NULL, referenceTypes = NULL) {
  VALID_MERGE_METHODS <- c("mean", "max")
  VALID_NORM_METHODS <- c("rank", "median", "mean")
  VALID_ENRICH_METHODS <- c("ORA", "GSEA")

  # Null Check
  analyteLists <- testNull(analyteLists)
  analyteListFiles <- testNull(analyteListFiles)
  analyteTypes <- testNull(analyteTypes)
  enrichMethod <- testNull(enrichMethod)
  organism <- testNull(organism)
  enrichDatabase <- testNull(enrichDatabase)
  enrichDatabaseFile <- testNull(enrichDatabaseFile)
  enrichDatabaseType <- testNull(enrichDatabaseType)
  enrichDatabaseDescriptionFile <- testNull(enrichDatabaseDescriptionFile)
  referenceLists <- testNull(referenceLists)
  referenceTypes <- testNull(referenceTypes)
  error_msg <- parameterErrorMessage(
    enrichMethod = enrichMethod, organism = organism, collapseMethod = collapseMethod, minNum = minNum, maxNum = maxNum,
    fdrMethod = fdrMethod, sigMethod = sigMethod, fdrThr = fdrThr, topThr = topThr, reportNum = reportNum, isOutput = isOutput,
    outputDirectory = outputDirectory, dagColor = dagColor, hostName = hostName, cache = cache
  )
  if (!is.null(error_msg)) {
    stop(error_msg)
  }

  # Verify parameters
  mergeMethod <- tolower(mergeMethod)
  normalizationMethod <- tolower(normalizationMethod)
  enrichMethod <- toupper(enrichMethod)
  if (enrichMethod == "ORA" && !isMetaAnalysis) {
    stop("ORA only supports meta-analysis. isMetaAnalysis must be set to TRUE")
  }
  if (!(mergeMethod %in% VALID_MERGE_METHODS)) {
    stop(paste0(mergeMethod, " is not a valid merge method.\nValid options are: ", paste(VALID_MERGE_METHODS, collapse = ", ")))
  }
  if (!(normalizationMethod %in% VALID_NORM_METHODS)) {
    stop(paste0(normalizationMethod, " is not a valid normalization method.\nValid options are: ", paste(VALID_NORM_METHODS, collapse = ", ")))
  }
  if (length(analyteLists) != length(analyteTypes) && length(analyteListFiles) != length(analyteTypes)) {
    stop("analyte lists and analyteTypes must be the same length.")
  }
  if (!(enrichMethod %in% VALID_ENRICH_METHODS)) {
    stop(paste0(enrichMethod, " is not a valid enrichment method for multiomics.\nValid options are: ", paste(VALID_ENRICH_METHODS, collapse = ", ")))
  }
  if (length(enrichDatabase) > 1 || length(enrichDatabaseFile) > 1) {
    stop("Only one enrichDatabase or enrichDatabaseFile can be specified for multiomics.")
  }

  if (enrichMethod == "ORA") {
    cat("Performing multi-omics ORA\nLoading the functional categories...\n")
    enrichD <- loadGeneSet(
      organism = organism, enrichDatabase = enrichDatabase, enrichDatabaseFile = enrichDatabaseFile, enrichDatabaseType = enrichDatabaseType,
      enrichDatabaseDescriptionFile = enrichDatabaseDescriptionFile, cache = cache, hostName = hostName
    )

    geneSet <- enrichD$geneSet
    geneSetDes <- enrichD$geneSetDes
    geneSetDag <- enrichD$geneSetDag
    geneSetNet <- enrichD$geneSetNet
    databaseStandardId <- enrichD$standardId
    rm(enrichD)

    cat("Loading the ID lists...\n")
    interest_lists <- list()
    if (is.null(analyteLists)) {
      for (i in analyteLists) {
        interestingGeneMap <- loadInterestGene(
          organism = organism, dataType = "list", inputGeneFile = analyteListFiles[i], inputGene = NULL,
          geneType = analyteTypes[i], collapseMethod = collapseMethod, cache = cache,
          hostName = hostName, geneSet = geneSet
        )
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

    ## Meta-analysis
  }
}
