#' @title Multi-omics ORA
#' @inheritParams WebGestaltRMultiOmics
WebGestaltRMultiOmicsOra <- function(analyteLists = NULL, analyteListFiles = NULL, analyteTypes = NULL, enrichMethod = "ORA", organism = "hsapiens",
                                     enrichDatabase = NULL, enrichDatabaseFile = NULL, enrichDatabaseType = NULL, enrichDatabaseDescriptionFile = NULL,
                                     collapseMethod = "mean", minNum = 10, maxNum = 500, fdrMethod = "BH", sigMethod = "fdr", fdrThr = 0.05,
                                     topThr = 10, reportNum = 100, setCoverNum = 10, perNum = 1000, gseaP = 1, isOutput = TRUE, outputDirectory = getwd(),
                                     projectName = NULL, dagColor = "binary", nThreads = 1, cache = NULL, hostName = "https://www.webgestalt.org/",
                                     useWeightedSetCover = TRUE, useAffinityPropagation = FALSE, usekMedoid = FALSE, kMedoid_k = 25,
                                     referenceLists = NULL, referenceListFiles = NULL, referenceTypes = NULL, listNames = null) {
  projectDir <- file.path(outputDirectory, paste0("Project_", projectName))
  cat("Performing multi-omics ORA\nLoading the functional categories...\n")
  all_sets <- .load_meta_gmt(enrichDatabase, enrichDatabaseFile, enrichDatabaseDescriptionFile, enrichDatabaseType, analyteLists, analyteListFiles, analyteTypes, organism, cache, hostName)
  cat("Loading the ID lists...\n")
  interest_lists <- list()
  interestGeneMaps  <- list()
  if (is.null(analyteLists)) {
    for (i in seq_along(analyteListFiles)) {
      interestingGeneMap <- loadInterestGene(
        organism = organism, dataType = "list", inputGeneFile = analyteListFiles[i], inputGene = NULL,
        geneType = analyteTypes[i], collapseMethod = collapseMethod, cache = cache,
        hostName = hostName, geneSet = all_sets[["geneSet"]][[i]]
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
        organism = organism, dataType = "list", inputGeneFile = NULL, inputGene = analyteLists[i],
        geneType = analyteTypes[i], collapseMethod = collapseMethod, cache = cache,
        hostName = hostName, geneSet = all_sets[["geneSet"]][[i]]
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

  # Load Gene Sets
  cat("Loading the reference lists...\n")
  reference_lists <- list()
  if (is.null(referenceLists)) {
    for (i in seq_along(referenceListFiles)) {
      referenceGeneList <- loadReferenceGene(
        organism = organism, referenceGeneFile = referenceListFiles[i],
        referenceGene = NULL, referenceGeneType = referenceTypes[i],
        referenceSet = NULL, collapseMethod = collapseMethod,
        hostName = hostName, geneSet = all_sets[["geneSet"]][[i]],
        interestGeneList = interest_lists[[i]],
        cache = cache
      )
      reference_lists[[i]] <- referenceGeneList
    }
  } else {
    for (i in seq_along(analyteLists)) {
      referenceGeneList <- loadReferenceGene(
        organism = organism, referenceGeneFile = NULL,
        referenceGene = referenceLists[i], referenceGeneType = NULL,
        referenceSet = NULL, collapseMethod = collapseMethod,
        hostName = hostName, geneSet = all_sets[["geneSet"]][[i]],
        interestGeneList = interest_lists[[i]],
        cache = cache
      )
      reference_lists[[i]] <- referenceGeneList
    }
  }
  cat("Running multi-omics ORA...\n")
  oraRes <- multiOraEnrichment(interest_lists, reference_lists, all_sets[["geneSet"]],
    minNum = minNum, maxNum = maxNum, fdrMethod = fdrMethod, sigMethod = sigMethod,
    fdrThr = fdrThr, topThr = topThr
  )
  if (is.null(oraRes)) {
    return(NULL)
  }
  cat("Generating the report...\n")
  dir.create(projectDir, showWarnings = FALSE)
  for (i in 2:length(oraRes$enriched)) {
    interestingGeneMap <- interestGeneMaps[[i - 1]]
    enrichedSig <- oraRes$enriched[[i]]
    insig <- oraRes$background[[i]]
    # geneSetDag <- all_sets[["geneSetDag"]][[i - 1]]
    geneSetDes <- all_sets[["geneSetDes"]][[i - 1]]
    geneSet <- all_sets[["geneSet"]][[i - 1]]
    clusters <- list()
    geneTables <- list()

    if (!is.null(enrichedSig)) {
      if (!is.null(geneSetDes)) { ####### Add extra description information ###########
        enrichedSig <- enrichedSig %>%
          left_join(geneSetDes, by = "geneSet") %>%
          select(.data$geneSet, .data$description, .data$link, .data$size, .data$overlap, .data$expect, .data$enrichmentRatio, .data$pValue, .data$FDR, .data$overlapId) %>%
          arrange(.data$FDR, .data$pValue, desc(.data$size)) %>%
          mutate(description = ifelse(is.na(.data$description), "", .data$description)) # now des could be mixture
      } else {
        enrichedSig <- enrichedSig %>%
          select(.data$geneSet, .data$link, .data$size, .data$overlap, .data$expect, .data$enrichmentRatio, .data$pValue, .data$FDR, .data$overlapId) %>%
          arrange(.data$FDR, .data$pValue, desc(.data$size))
      }

      geneTables <- getGeneTables(organism, enrichedSig, "overlapId", interestingGeneMap)
      if (organism != "others") {
        enrichedSig$link <- mapply(
          function(link, geneList) linkModification("ORA", link, geneList, interestingGeneMap, hostName),
          enrichedSig$link,
          enrichedSig$overlapId
        )
      }

      if ("database" %in% colnames(geneSet)) {
        # Add source database for multiple databases
        enrichedSig <- enrichedSig %>% left_join(unique(geneSet[, c("geneSet", "database")]), by = "geneSet")
      }
      if (organism != "others" && analyteTypes[[i - 1]] != interestStandardId) {
        outputEnrichedSig <- mapUserId(enrichedSig, "overlapId", interestingGeneMap)
      } else {
        outputEnrichedSig <- enrichedSig
      }

      if (isOutput) {
        write_tsv(outputEnrichedSig, file.path(projectDir, paste0("enrichment_results_", projectName, "_", listNames[i - 1], ".txt")))
        idsInSet <- sapply(enrichedSig$overlapId, strsplit, split = ";")
        names(idsInSet) <- enrichedSig$geneSet
        minusLogP <- -log(enrichedSig$pValue)
        minusLogP[minusLogP == Inf] <- -log(.Machine$double.eps)
        apRes <- NULL
        wscRes <- NULL
        kRes <- NULL
        if (useAffinityPropagation) {
          apRes <- affinityPropagation(idsInSet, minusLogP)
        }
        if (useWeightedSetCover) {
          wscRes <- weightedSetCover(idsInSet, 1 / minusLogP, setCoverNum, nThreads)
        }
        if (usekMedoid) {
          kRes <- kMedoid(idsInSet, minusLogP, maxK = kMedoid_k)
        }
        if (!is.null(apRes)) {
          writeLines(sapply(apRes$clusters, paste, collapse = "\t"), file.path(projectDir, paste0("enriched_geneset_ap_clusters_", projectName, "_", listNames[i - 1], ".txt")))
        } else {
          apRes <- NULL
        }
        clusters$ap <- apRes
        if (!is.null(kRes)) {
          writeLines(sapply(kRes$clusters, paste, collapse = "\t"), file.path(projectDir, paste0("enriched_geneset_kmedoid_clusters_", projectName, "_", listNames[i - 1], ".txt")))
        } else {
          kRes <- NULL
        }
        clusters$km <- kRes
        if (!is.null(wscRes$topSets)) {
          writeLines(c(paste0("# Coverage: ", wscRes$coverage), wscRes$topSets), file.path(projectDir, paste0("enriched_geneset_wsc_topsets_", projectName, "_", listNames[i - 1], ".txt")))
          clusters$wsc <- list(representatives = wscRes$topSets, coverage = wscRes$coverage)
        } else {
          clusters$wsc <- NULL
        }
      }
    }
  }
}
