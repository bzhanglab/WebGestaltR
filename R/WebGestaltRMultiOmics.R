#' @title WebGestaltRMultiOmics
#' @description Perform multi-omics analysis using WebGestaltR
#' @param analyteLists A list of analyte lists
#' @param analyteListFiles A list of analyte list files
#' @param analyteTypes A list of analyte types
#' @param enrichMethod Enrichment method, either \code{ORA} or \code{GSEA}
#' @param organism The organism to use
#' @param enrichDatabase The database to use
#' @param enrichDatabaseFile The database file to use
#' @param enrichDatabaseType The database type to use
#' @param enrichDatabaseDescriptionFile The database description file to use
#' @export
WebGestaltRMultiOmics <- function(analyteLists = NULL, analyteListFiles = NULL, analyteTypes = NULL, enrichMethod = "ORA", organism = "hsapiens", enrichDatabase = NULL,
                                  enrichDatabaseFile = NULL, enrichDatabaseType = NULL, enrichDatabaseDescriptionFile = NULL,
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
  referenceLists <- testNull(referenceLists)
  referenceTypes <- testNull(referenceTypes)
  enrichMethod <- testNull(enrichMethod)
  organism <- testNull(organism)
  enrichDatabase <- testNull(enrichDatabase)
  enrichDatabaseFile <- testNull(enrichDatabaseFile)
  enrichDatabaseType <- testNull(enrichDatabaseType)
  enrichDatabaseDescriptionFile <- testNull(enrichDatabaseDescriptionFile)
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
    databases <- c()
    if (!is.null(enrichDatabase)) { # Need to get correct name for metabolite databases
      if (length(unique(analyteTypes)) == 1) {
        databases <- enrichDatabase
      } else {
        types_processed <- c()
        for (i in seq_along(analyteTypes)) {
          if (analyteTypes[i] %in% types_processed) {
            next
          }
          databases <- c(databases, get_gmt_file(hostName, analyteTypes[i], enrichDatabase[i], organism, cache))
          types_processed <- c(types_processed, analyteTypes[i])
        }
        databases <- unique(databases)
      }
    } else {
      databases <- NULL
    }
    all_sets <- loadGeneSet(
      organism = organism, enrichDatabase = databases, enrichDatabaseFile = enrichDatabaseFile, enrichDatabaseType = enrichDatabaseType,
      enrichDatabaseDescriptionFile = enrichDatabaseDescriptionFile, cache = cache, hostName = hostName
    )
    if (length(all_sets) > 1) {
      geneSet <- all_sets[[1]]$geneSet
      geneSetDes <- all_sets[[1]]$geneSetDes
      geneSetNet <- all_sets[[1]]$geneSetNet
      for (i in 2:length(all_sets)) {
        geneSet <- rbind(geneSet, all_sets[[i]]$geneSet)
        geneSetDes <- rbind(geneSetDes, all_sets[[i]]$geneSetDes)
        geneSetDag <- rbind(geneSetDag, all_sets[[i]]$geneSetDag)
        geneSetNet <- rbind(geneSetNet, all_sets[[i]]$geneSetNet)
        databaseStandardId <- "multiomics"
      }
    } else {
      geneSet <- all_sets$geneSet
      geneSetDag <- all_sets$geneSetDag
      geneSetNet <- all_sets$geneSetNet
      databaseStandardId <- all_sets$standardId
    }

    rm(all_sets)

    cat("Loading the ID lists...\n")
    interest_lists <- list()
    if (is.null(analyteLists)) {
      for (i in seq_along(analyteListFiles)) {
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
    } else {
      for (i in seq_along(analyteLists)) {
        interestingGeneMap <- loadInterestGene(
          organism = organism, dataType = "list", inputGeneFile = analyteLists[i], inputGene = NULL,
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

    # Load Gene Sets
    cat("Loading the reference lists...\n")

    ## Meta-analysis
  }
}
