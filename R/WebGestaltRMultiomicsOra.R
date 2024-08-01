#' @title Multi-omics ORA
#' importFrom dplyr bind_rows left_join arrange select desc %>% mutate distinct
#' importFrom readr write_tsv
#' @inheritParams WebGestaltRMultiOmics
WebGestaltRMultiOmicsOra <- function(analyteLists = NULL, analyteListFiles = NULL, analyteTypes = NULL, enrichMethod = "ORA", organism = "hsapiens",
                                     enrichDatabase = NULL, enrichDatabaseFile = NULL, enrichDatabaseType = NULL, enrichDatabaseDescriptionFile = NULL,
                                     collapseMethod = "mean", minNum = 10, maxNum = 500, fdrMethod = "BH", sigMethod = "fdr", fdrThr = 0.05,
                                     topThr = 10, reportNum = 100, setCoverNum = 10, perNum = 1000, gseaP = 1, isOutput = TRUE, outputDirectory = getwd(),
                                     projectName = NULL, dagColor = "binary", nThreads = 1, cache = NULL, hostName = "https://www.webgestalt.org/",
                                     useWeightedSetCover = TRUE, useAffinityPropagation = FALSE, usekMedoid = FALSE, kMedoid_k = 25,
                                     referenceLists = NULL, referenceListFiles = NULL, referenceTypes = NULL, referenceSets = NULL, listNames = NULL) {
  projectDir <- file.path(outputDirectory, paste0("Project_", projectName))
  cat("Performing multi-omics ORA\nLoading the functional categories...\n")
  all_sets <- .load_meta_gmt(enrichDatabase, enrichDatabaseFile, enrichDatabaseDescriptionFile, enrichDatabaseType, analyteLists, analyteListFiles, analyteTypes, organism, cache, hostName)
  cat("Loading the ID lists...\n")
  interest_lists <- list()
  interestGeneMaps <- list()
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

  cat("Loading the reference lists...\n")
  reference_lists <- list()
  if (!is.null(referenceListFiles)) {
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
  } else if (!is.null(referenceLists)) {
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
  } else { # use pre-defined reference lists
    for (i in seq_along(interest_lists)) {
      referenceGeneList <- loadReferenceGene(
        organism = organism, referenceGeneFile = NULL,
        referenceGene = NULL, referenceGeneType = NULL,
        referenceSet = referenceSets[[i]], collapseMethod = collapseMethod,
        hostName = hostName, geneSet = all_sets[["geneSet"]][[i]],
        interestGeneList = interest_lists[[i]],
        cache = cache
      )
      reference_lists[[i]] <- referenceGeneList
    }
  }
  cat("Running multi-omics ORA...\n")
  oraRes <- multiOraEnrichment(interest_lists, reference_lists, all_sets[["geneSet"]],
    minNum = minNum, maxNum = maxNum, fdrMethod = fdrMethod,
    sigMethod = sigMethod, fdrThr = fdrThr, topThr = topThr
  )
  if (is.null(oraRes)) {
    return(NULL)
  }
  cat("Generating the report...\n")
  geneSetDags <- all_sets[["geneSetDag"]]
  geneSetNets <- all_sets[["geneSetNet"]]
  enrichSigs <- list()
  enrichSigs[[1]] <- NULL
  clusters_list <- list()
  dir.create(projectDir, showWarnings = FALSE)
  insig_lists <- list()
  insig_lists[[1]] <- NULL
  geneTables_list <- list()
  enrichedSigs <- list()
  geneTables_list[[1]] <- NULL

  ## Meta-analysis

  for (i in 2:(length(oraRes$enriched) + 1)) {
    print(paste0("Processing list ", i))
    if (i == (length(oraRes$enriched) + 1)) {
      i <- 1
    }
    enrichedSig <- oraRes$enriched[[i]]
    if (i == 1) {
      interestingGeneMap <- list(mapped = interestGeneMaps[[1]]$mapped, unmapped = interestGeneMaps[[1]]$unmapped, standardId = interestGeneMaps[[1]]$standardId)
      for (j in seq_along(interestGeneMaps)) {
        if (j == 1) {
          next
        }
        old_names <- names(interestGeneMaps[[j]][["mapped"]])
        names(interestGeneMaps[[j]][["mapped"]]) <- names(interestGeneMaps[[1]][["mapped"]])
        interestingGeneMap[["mapped"]] <- rbind(interestingGeneMap[["mapped"]], interestGeneMaps[[j]][["mapped"]])
        interestingGeneMap[["unmapped"]] <- append(unlist(interestingGeneMap[["unmapped"]]), unlist(interestGeneMaps[[j]][["unmapped"]]))
        names(interestGeneMaps[[j]][["mapped"]]) <- old_names
      }
      if ("geneSetDes" %in% names(all_sets)) {
        geneSetDes <- all_sets[["geneSetDes"]][[1]]
        geneSet <- all_sets[["geneSet"]][[1]]
      } else {
        geneSetDes <- NULL
        geneSet <- all_sets[["geneSet"]][[1]]
      }
      for (j in seq_along(all_sets[["geneSet"]])) {
        if (j == 1) {
          next
        }
        geneSet <- rbind(geneSet, all_sets[["geneSet"]][[j]])
        if (!is.null(geneSetDes)) {
          if (length(all_sets[["geneSetDes"]]) >= (j) && !is.null(all_sets[["geneSetDes"]][[j]]) && (ncol(all_sets[["geneSetDes"]][[j]]) == ncol(geneSetDes))) {
            geneSetDes <- rbind(geneSetDes, all_sets[["geneSetDes"]][[j]])
          }
        }
      }
      enrichedSig <- enrichedSig %>%
        left_join(geneSet[, c("geneSet", "description")], by = "geneSet")
      names(enrichedSig)[names(enrichedSig) == "description"] <- "link"
    } else {
      interestingGeneMap <- interestGeneMaps[[i - 1]]
      if ("geneSetDes" %in% names(all_sets)) {
        if (length(all_sets[["geneSetDes"]]) >= (i - 1)) {
          geneSetDes <- all_sets[["geneSetDes"]][[i - 1]]
          geneSet <- all_sets[["geneSet"]][[i - 1]]
        } else {
          geneSetDes <- NULL
          geneSet <- all_sets[["geneSet"]][[i - 1]]
        }
      } else {
        geneSetDes <- NULL
        geneSet <- all_sets[["geneSet"]][[i - 1]]
      }
    }
    enrichedSigs[[i]] <- enrichedSig
    insig <- oraRes$background[[i]]
    insig_lists[[i]] <- insig
    clusters <- list()
    if (!is.null(enrichedSig)) {
      if (!is.null(geneSetDes)) { ####### Add extra description information ###########
        enrichedSig <- enrichedSig %>%
          left_join(geneSetDes, by = "geneSet") %>%
          select(.data$geneSet, .data$description, .data$link, .data$size, .data$overlap, .data$expect, .data$enrichmentRatio, .data$pValue, .data$FDR, .data$overlapId) %>%
          arrange(.data$FDR, .data$pValue, desc(.data$size)) %>%
          mutate(description = ifelse(is.na(.data$description), "", .data$description)) %>%
          distinct(.data$geneSet, .keep_all = TRUE)
      } else {
        enrichedSig <- enrichedSig %>%
          select(.data$geneSet, .data$link, .data$size, .data$overlap, .data$expect, .data$enrichmentRatio, .data$pValue, .data$FDR, .data$overlapId) %>%
          arrange(.data$FDR, .data$pValue, desc(.data$size)) %>%
          distinct(.data$geneSet, .keep_all = TRUE)
      }
      if (i == 1) {
        print(enrichedSig)
      }
      geneTables <- getGeneTables(organism, enrichedSig, "overlapId", interestingGeneMap)
      geneTables_list[[i]] <- geneTables
      if (organism != "others" && i != 1) {
        genelist_copy <- enrichedSig$overlapId
        enrichedSig$link <- mapply(
          function(link, geneList) linkModification("ORA", link, geneList, interestingGeneMap, hostName),
          enrichedSig$link,
          genelist_copy
        )
      } else if (organism != "others") {
        idsInSet <- sapply(enrichedSig$overlapId, strsplit, split = ";")
        names(idsInSet) <- enrichedSig$geneSet
        old_links <- enrichedSig$link
        for (k in seq_along(enrichedSig$link)) {
          new_link <- metaLinkModification("ORA", enrichedSig$link[[k]], unlist(idsInSet[[k]]), interestGeneMaps, hostName, enrichedSig$geneSet[[k]])
          if (!is.null(new_link)) {
            enrichedSig$link[[k]] <- new_link
          } else {
            enrichedSig$link[[k]] <- old_links[[k]]
          }
        }
      }

      if ("database" %in% colnames(geneSet)) {
        # Add source database for multiple databases
        enrichedSig <- enrichedSig %>% left_join(unique(geneSet[, c("geneSet", "database")]), by = "geneSet")
      }
      if (i == 1) {
        outputEnrichedSig <- enrichedSig
      } else if (organism != "others" && analyteTypes[[i - 1]] != interestStandardId) {
        outputEnrichedSig <- mapUserId(enrichedSig, "overlapId", interestingGeneMap)
      } else {
        outputEnrichedSig <- enrichedSig
      }


      if (isOutput) {
        if (i == 1) {
          write_tsv(outputEnrichedSig, file.path(projectDir, paste0("enrichment_results_", projectName, ".txt")))
        } else {
          write_tsv(outputEnrichedSig, file.path(projectDir, paste0("enrichment_results_", projectName, "_", listNames[i - 1], ".txt")))
        }
        idsInSet <- sapply(enrichedSig$overlapId, strsplit, split = ";")
        names(idsInSet) <- enrichedSig$geneSet
        minusLogP <- -log(enrichedSig$pValue)
        minusLogP[minusLogP == Inf] <- -log(.Machine$double.eps)
        apRes <- NULL
        wscRes <- NULL
        kRes <- NULL
        if (i == 1) {
          name_ending <- projectName
        } else {
          name_ending <- paste0(projectName, "_", listNames[i - 1])
        }
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
          writeLines(sapply(apRes$clusters, paste, collapse = "\t"), file.path(projectDir, paste0("enriched_geneset_ap_clusters_", name_ending, ".txt")))
        } else {
          apRes <- NULL
        }
        clusters$ap <- apRes
        if (!is.null(kRes)) {
          writeLines(sapply(kRes$clusters, paste, collapse = "\t"), file.path(projectDir, paste0("enriched_geneset_kmedoid_clusters_", name_ending, ".txt")))
        } else {
          kRes <- NULL
        }
        clusters$km <- kRes
        if (!is.null(wscRes$topSets)) {
          writeLines(c(paste0("# Coverage: ", wscRes$coverage), wscRes$topSets), file.path(projectDir, paste0("enriched_geneset_wsc_topsets_", name_ending, ".txt")))
          clusters$wsc <- list(representatives = wscRes$topSets, coverage = wscRes$coverage)
        } else {
          clusters$wsc <- NULL
        }
        clusters_list[[i]] <- clusters
      }
    }
    if (i == 1) {
      enrichedSig <- enrichedSig %>% distinct(.data$geneSet, .keep_all = TRUE)
    }
    enrichSigs[[i]] <- enrichedSig
    if (i == 1) {
      break
    }
  }
  if (isOutput) {
    ############## Create report ##################
    cat("Generate the final report...\n")
    createMetaReport(
      hostName = hostName, outputDirectory = outputDirectory, organism = organism,
      projectName = projectName, enrichMethod = enrichMethod, geneSet_list = all_sets[["geneSet"]],
      geneSetDes_list = geneSetDes, geneSetDag_list = geneSetDags, geneSetNet_list = geneSetNets,
      interestingGeneMap_list = interestGeneMaps, referenceGeneList_list = reference_lists,
      enrichedSig_list = enrichSigs, background_list = insig_lists, geneTables_list = geneTables_list,
      clusters_list = clusters_list, enrichDatabase_list = enrichDatabase,
      enrichDatabaseFile_list = enrichDatabaseFile, enrichDatabaseType_list = enrichDatabaseType,
      enrichDatabaseDescriptionFile_list = enrichDatabaseDescriptionFile,
      interestGeneFile_list = analyteListFiles, interestGene_list = interest_lists,
      interestGeneType_list = analyteTypes, collapseMethod = collapseMethod,
      referenceGeneFile_list = referenceListFiles, referenceGene_list = referenceLists,
      referenceGeneType_list = referenceTypes, referenceSet_list = referenceLists, minNum = minNum,
      maxNum = maxNum, fdrMethod = fdrMethod, sigMethod = sigMethod, fdrThr = fdrThr,
      topThr = topThr, reportNum = reportNum, dagColor = dagColor, listNames = listNames
    )

    cwd <- getwd()
    setwd(projectDir)
    zip(paste0("Project_", projectName, ".zip"), ".", flags = "-rq")
    setwd(cwd)

    cat("Results can be found in the ", projectDir, "!\n", sep = "")
  }
  return(outputEnrichedSig)
}
