mapUserId <- function(enrichedSig, geneColumn, interestingGeneMap) {
    #### map entrez gene back to the original user id and add one more column to the enrichedSig
    standardId <- interestingGeneMap$standardId

    mapgene <- interestingGeneMap$mapped[, c("userId", standardId)]
    gene <- enrichedSig[[geneColumn]]
    gene <- strsplit(gene, ";")
    gene <- unlist(lapply(gene, geneM, mapgene))
    # print(gene)
    # if (standardId == "rampc"){
    # 	gene <- unlist(lapply(gene, function(x) {
    # 		gene_array <- strsplit(x,";")
    # 		if (is.na(x)) {
    # 			return(NA);
    # 		} else {
    # 			return(paste(strsplit(x, ":")[-1], sep=":"))
    # 		}
    # 	}))
    # }
    enrichedSig <- data.frame(enrichedSig, userId = gene, stringsAsFactors = FALSE)
    return(enrichedSig)
}

geneM <- function(geneList, mappingTable) {
    if (length(geneList) == 1 && is.na(geneList)) {
        ### The categories outputted from GSEA may not have leading edge genes. TODO: obsolete
        return(NA)
    } else {
        u <- mappingTable[mappingTable[[2]] %in% geneList, ][[1]]
        # although user ID could contain ;, like in some gene symbols.
        # but this is only concatenated in output
        u <- paste(u, collapse = ";")
        return(u)
    }
}

#' @importFrom dplyr select
getMetaGSEAGeneTables <- function(organism, enrichedSig, interestingGeneMaps, listNames) {
    table <- list()
    for (i in seq_len(nrow(enrichedSig))) {
        genes <- unlist(unique(unlist(strsplit(enrichedSig[[i, "leadingEdgeId"]], ";"))))
        geneSetId <- enrichedSig[[i, "geneSet"]]
        if (length(genes) == 1 && is.na(genes)) {
            table[[geneSetId]] <- list()
        } else {
            if (organism != "others") {
                row <- data.frame("Analyte" = genes)
                new_genes <- list()
                score_lists <- list()
                for (j in seq_along(genes)) {
                    gene <- genes[j]
                    new_label <- ""
                    for (k in seq_along(interestingGeneMaps)) {
                        mapping <- select(interestingGeneMaps[[k]]$mapped, .data$geneSymbol, .data$score, interestingGeneMaps[[k]]$standardId)

                        if (gene %in% mapping[[interestingGeneMaps[[k]]$standardId]]) {
                            score_lists[[k]] <- unlist(mapping[mapping[[interestingGeneMaps[[k]]$standardId]] %in% gene, "score"])
                            new_label <- unlist(mapping[mapping[[interestingGeneMaps[[k]]$standardId]] %in% gene, "geneSymbol"])[1]
                        } else {
                            score_lists[[k]] <- "-"
                        }
                    }
                    if (new_label == "") {
                        new_label <- gene
                    }
                    new_genes[[j]] <- new_label
                }
                row <- data.frame("Analyte" = new_genes)
                for (k in seq_along(interestingGeneMaps)) {
                    row[[listNames[[k]]]] <- score_lists[[k]]
                }
                table[[geneSetId]] <- row
            } else {
                table[[geneSetId]] <- data.frame("Analyte" = genes)
            }
        }
    }
    return(table)
}

#' @importFrom dplyr select
getGeneTables <- function(organism, enrichedSig, geneColumn, interestingGeneMap) {
    if (organism != "others") {
        standardId <- interestingGeneMap$standardId
        mapping <- select(interestingGeneMap$mapped, .data$userId, .data$geneSymbol, .data$geneName, .data$gLink, standardId)
        if ("score" %in% colnames(interestingGeneMap$mapped)) {
            mapping$score <- interestingGeneMap$mapped$score
        }
    }
    table <- list()
    for (i in seq_len(nrow(enrichedSig))) {
        genes <- enrichedSig[[i, geneColumn]]
        geneSetId <- enrichedSig[[i, "geneSet"]]
        if (length(genes) == 1 && is.na(genes)) {
            table[[geneSetId]] <- list()
        } else {
            genes <- unlist(strsplit(genes, ";"))
            if (organism != "others") {
                table[[geneSetId]] <- mapping[mapping[[standardId]] %in% genes, ]
            } else {
                table[[geneSetId]] <- data.frame("userId" = genes)
            }
        }
    }
    return(table)
}

#' @importFrom dplyr filter bind_rows
getTopGseaResults <- function(results, topThr) {
    if (is.wholenumber(topThr)) {
        posThr <- topThr
        negThr <- topThr
    } else {
        posThr <- floor(topThr) + 1
        negThr <- floor(topThr)
    }
    posRes <- filter(results, .data$normalizedEnrichmentScore > 0)
    if (nrow(posRes) > posThr) {
        posSig <- posRes[1:posThr, ]
        posInsig <- posRes[(posThr + 1):nrow(posRes), ]
    } else {
        posSig <- posRes
        posInsig <- NULL
    }
    negRes <- filter(results, .data$normalizedEnrichmentScore < 0)
    if (nrow(negRes) > negThr) {
        negSig <- negRes[1:negThr, ]
        negInsig <- negRes[(negThr + 1):nrow(negRes), ]
    } else {
        negSig <- negRes
        negInsig <- NULL
    }
    sig <- bind_rows(posSig, negSig)
    insig <- bind_rows(posInsig, negInsig)
    if (nrow(sig) == 0) {
        sig <- NULL
    }
    if (nrow(insig) == 0) {
        insig <- NULL
    }
    return(list(sig, insig))
}

#' @importFrom dplyr filter bind_rows
getTopMetaGseaResults <- function(results, topThr) {
    if (is.wholenumber(topThr)) {
        topThr <- topThr
    } else {
        topThr <- floor(topThr) + 1
    }
    if (nrow(results) > topThr) {
        sig <- results[1:topThr, ]
        insig <- results[(topThr + 1):nrow(results), ]
    } else {
        sig <- results
        insig <- NULL
    }
    return(list(sig, insig))
}

#' keepRep
#'
#' Add representatives of redundancy-reduced clusters to
#' topResult if they are missing.
#'
#' @keywords internal
#'
keepRep <- function(topResult, allResult, reps) {
    missing <- NULL
    for (rep in reps) {
        if (!rep %in% topResult$geneSet) {
            missing <- c(missing, rep)
        }
    }
    if (!is.null(missing)) {
        topResult <- rbind(topResult, allResult[allResult$geneSet %in% missing, ])
    }
    return(topResult)
}

removeFileProtocol <- function(path) {
    return(normalizePath(sub("^file://", "", path), mustWork = FALSE))
}

sanitizeFileName <- function(name) {
    return(gsub("[[:punct:]]", "_", name))
}
