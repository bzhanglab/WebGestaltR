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
    all_sets <- .load_meta_gmt(
        enrichDatabase, enrichDatabaseFile, enrichDatabaseDescriptionFile, enrichDatabaseType, analyteLists, analyteListFiles, analyteTypes,
        organism, cache, hostName
    )
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
                interestGeneList <- interestingGeneMap$mapped %>%
                    select(interestStandardId, .data$score) %>%
                    distinct()
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
                interestGeneList <- interestingGeneMap$mapped %>%
                    select(interestStandardId, .data$score) %>%
                    distinct()
            }
            interest_lists[[i]] <- interestGeneList
        }
    }
    cat("Running multi-omics GSEA...\n")

    gseaRes <- multiGseaEnrichment(
        hostName = hostName, outputDirectory = outputDirectory, projectName = projectName, geneRankList_list = interest_lists, geneSet_list = all_sets[["geneSet"]],
        geneSetDes_list = all_sets[["geneSetDes"]], collapseMethod = "mean", minNum = minNum, maxNum = maxNum, sigMethod = sigMethod, fdrThr = fdrThr,
        topThr = topThr, perNum = perNum, p = gseaP, isOutput = isOutput, saveRawGseaResult = saveRawGseaResult, plotFormat = gseaPlotFormat,
        nThreads = nThreads, listNames = listNames
    )

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
    geneTables_list[[i]] <- NULL
    enrichedSig_list <- list()

    ## Meta-analysis

    for (i in 2:(length(gseaRes$enriched) + 1)) {
        print(paste0("Processing ", i, " list..."))
        if (i == (length(gseaRes$enriched) + 1)) {
            print("Meta-analysis")
            i <- 1
        }
        enrichedSig <- gseaRes$enriched[[i]]
        insig <- gseaRes$background[[i]]
        print("first")
        if (i == 1) {
            interestingGeneMap <- list(mapped = interestGeneMaps[[1]]$mapped, unmapped = interestGeneMaps[[1]]$unmapped, standardId = interestGeneMaps[[1]]$standardId)
            for (j in seq_along(interestGeneMaps)) {
                if (j == 1) {
                    next
                }
                old_names <- names(interestGeneMaps[[j]][["mapped"]])
                names(interestGeneMaps[[j]][["mapped"]]) <- names(interestGeneMaps[[1]][["mapped"]])
                interestingGeneMap[["mapped"]] <- rbind(interestingGeneMap[["mapped"]], interestGeneMaps[[j]][["mapped"]])
                interestingGeneMap[["unmapped"]] <- append(interestingGeneMap[["unmapped"]], interestGeneMaps[[j]][["unmapped"]])
                names(interestGeneMaps[[j]][["mapped"]]) <- old_names
            }
            if ("geneSetDes" %in% names(all_sets)) {
                if (length(all_sets[["geneSetDes"]]) < 1) {
                    geneSetDes <- all_sets[["geneSet"]][[1]]
                    geneSet <- all_sets[["geneSet"]][[1]]
                } else {
                    geneSetDes <- all_sets[["geneSetDes"]][[1]]
                    geneSet <- all_sets[["geneSet"]][[1]]
                    for (j in seq_along(all_sets[["geneSet"]])) {
                        if (j == 1) {
                            next
                        }
                        geneSet <- rbind(geneSet, all_sets[["geneSet"]][[j]])
                        if (length(all_sets[["geneSetDes"]]) >= (j) && !is.null(all_sets[["geneSetDes"]][[j]]) && (ncol(all_sets[["geneSetDes"]][[j]]) == ncol(geneSetDes))) {
                            geneSetDes <- rbind(geneSetDes, all_sets[["geneSetDes"]][[j]])
                        }
                    }
                }
            } else {
                geneSetDes <- all_sets[["geneSet"]][[1]]
                geneSet <- all_sets[["geneSet"]][[1]]
                for (j in seq_along(all_sets[["geneSet"]])) {
                    if (j == 1) {
                        next
                    }
                    geneSet <- rbind(geneSet, all_sets[["geneSet"]][[j]])
                }
            }
            geneSet <- geneSet %>% distinct(.data$geneSet, .keep_all = TRUE)
            geneSetDes <- geneSetDes %>% distinct(.data$geneSet, .keep_all = TRUE)
        } else {
            interestingGeneMap <- interestGeneMaps[[i - 1]]
            if ("geneSetDes" %in% names(all_sets)) {
                if (length(all_sets[["geneSetDes"]]) >= (i - 1) && !is.null(all_sets[["geneSetDes"]][[i - 1]])) {
                    geneSetDes <- all_sets[["geneSetDes"]][[i - 1]]
                    geneSet <- all_sets[["geneSet"]][[i - 1]]
                } else {
                    geneSet <- NULL
                    geneSet <- all_sets[["geneSet"]][[i - 1]]
                }
            } else {
                geneSet <- NULL
                geneSet <- all_sets[["geneSet"]][[i - 1]]
            }
        }
        insig_lists[[i]] <- insig

        clusters <- list()
        geneTables <- list()
        if (!is.null(enrichedSig)) {
            if (!is.null(geneSetDes)) { ####### Add extra description information ###########
                print("here")
                enrichedSig <- enrichedSig %>%
                    left_join(geneSetDes, by = "geneSet") %>%
                    select(.data$geneSet, .data$description, .data$link, .data$enrichmentScore, .data$normalizedEnrichmentScore, .data$pValue, .data$FDR, .data$size, .data$plotPath, .data$leadingEdgeNum, .data$leadingEdgeId) %>%
                    arrange(.data$FDR, .data$pValue, desc(.data$normalizedEnrichmentScore)) %>%
                    mutate(description = ifelse(is.na(.data$description), "", .data$description))
            } else {
                enrichedSig <- enrichedSig %>%
                    select(.data$geneSet, .data$link, .data$enrichmentScore, .data$normalizedEnrichmentScore, .data$pValue, .data$FDR, .data$size, .data$plotPath, .data$leadingEdgeNum, .data$leadingEdgeId) %>%
                    arrange(.data$FDR, .data$pValue, desc(.data$normalizedEnrichmentScore))
            }

            enrichedSig <- enrichedSig %>% distinct(.data$geneSet, .keep_all = TRUE)
            if (i != 1) {
                enrichedSig_list[[i - 1]] <- enrichedSig
            }
            print("next")

            if (organism != "others" && i != 1) {
                geneTables <- getGeneTables(organism, enrichedSig, "leadingEdgeId", interestingGeneMap)
                geneTables_list[[i]] <- geneTables
                enrichedSig$link <- mapply(
                    function(link, geneList) linkModification("GSEA", link, geneList, interestingGeneMap, hostName),
                    enrichedSig$link,
                    enrichedSig$leadingEdgeId
                )
            } else if (organism != "others") {
                geneTables <- getGeneTables(organism, enrichedSig, "leadingEdgeId", interestingGeneMap)
                geneTables_list[[i]] <- geneTables
                metaGeneTables <- getMetaGSEAGeneTables(organism, enrichedSig, interestGeneMaps, listNames)
                for (k in seq_along(enrichedSig$link)) {
                    old_link <- enrichedSig$link[[k]]
                    tryCatch(
                        {
                            new_link <- metaLinkModification("GSEA", enrichedSig$link[[k]], strsplit(enrichedSig$leadingEdgeId[[k]], ";"), interestGeneMaps, hostName, enrichedSig$geneSet[[k]])
                            if (!is.null(new_link) && !is.na(new_link) && new_link != "") {
                                enrichedSig$link[[k]] <- new_link
                            } else {
                                enrichedSig$link[[k]] <- old_link
                            }
                        },
                        error = function(e) {
                            enrichedSig$link[[k]] <- old_link
                        },
                        warning = function(w) {
                            enrichedSig$link[[k]] <- old_link
                        }
                    )
                }
            }
            if ("database" %in% colnames(geneSet)) {
                # add source database for multiple databases
                enrichedSig <- enrichedSig %>% left_join(unique(geneSet[, c("geneSet", "database")]), by = "geneSet")
            }
            if (i == 1) {
                outputEnrichedSig <- enrichedSig
            } else if (organism != "others" && analyteTypes[[i - 1]] != interestStandardId) {
                outputEnrichedSig <- mapUserId(enrichedSig, "leadingEdgeId", interestingGeneMap)
            } else {
                outputEnrichedSig <- enrichedSig
            }
            if (isOutput) {
                write_tsv(outputEnrichedSig, file.path(projectDir, paste0("enrichment_results_", projectName, ".txt")))
                idsInSet <- sapply(enrichedSig$leadingEdgeId, strsplit, split = ";")
                names(idsInSet) <- enrichedSig$geneSet
                pValue <- enrichedSig$pValue
                pValue[pValue == 0] <- .Machine$double.eps
                signedLogP <- -log(pValue) * sign(enrichedSig$enrichmentScore)
                apRes <- NULL
                wscRes <- NULL
                kRes <- NULL
                if (useAffinityPropagation) {
                    apRes <- affinityPropagation(idsInSet, signedLogP)
                }
                if (useWeightedSetCover) {
                    wscRes <- weightedSetCover(idsInSet, 1 / signedLogP, setCoverNum, nThreads)
                }
                if (usekMedoid) {
                    kRes <- kMedoid(idsInSet, signedLogP, maxK = kMedoid_k)
                }
                if (!is.null(apRes)) {
                    tryCatch(
                        {
                            writeLines(sapply(apRes$clusters, paste, collapse = "\t"), file.path(projectDir, paste0("enriched_geneset_ap_clusters_", projectName, ".txt")))
                        },
                        error = function(e) {
                            cat("Error in writing ap clusters.\n")
                        }
                    )
                } else {
                    apRes <- NULL
                }
                clusters$ap <- apRes
                if (!is.null(kRes)) {
                    tryCatch(
                        {
                            writeLines(sapply(kRes$clusters, paste, collapse = "\t"), file.path(projectDir, paste0("enriched_geneset_kmedoid_clusters_", projectName, ".txt")))
                        },
                        error = function(e) {
                            cat("Error in writing kmedoid clusters.\n")
                        }
                    )
                } else {
                    kRes <- NULL
                }
                clusters$km <- kRes
                if (!is.null(wscRes$topSets)) {
                    writeLines(c(paste0("# Coverage: ", wscRes$coverage), wscRes$topSets), file.path(projectDir, paste0("enriched_geneset_wsc_topsets_", projectName, ".txt")))
                    clusters$wsc <- list(representatives = wscRes$topSets, coverage = wscRes$coverage)
                } else {
                    clusters$wsc <- NULL
                }
            }
        }
        clusters_list[[i]] <- clusters
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
            geneSetDes_list = all_sets[["geneSetDes"]], geneSetDag_list = geneSetDags, geneSetNet_list = geneSetNets,
            interestingGeneMap_list = interestGeneMaps, enrichedSig_list = enrichSigs, background_list = insig_lists,
            geneTables_list = geneTables_list, clusters_list = clusters_list, enrichDatabase_list = enrichDatabase,
            enrichDatabaseFile_list = enrichDatabaseFile, enrichDatabaseType_list = enrichDatabaseType,
            enrichDatabaseDescriptionFile_list = enrichDatabaseDescriptionFile,
            interestGeneFile_list = analyteListFiles, interestGene_list = interest_lists,
            interestGeneType_list = analyteTypes, collapseMethod = collapseMethod, minNum = minNum, perNum = perNum, p = gseaP,
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
