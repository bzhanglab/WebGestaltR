#' @importFrom dplyr filter select left_join mutate arrange %>% group_by inner_join
#' @importFrom stats p.adjust phyper
multiOraEnrichment <- function(interestGene, referenceGene, geneSet, minNum = 10,
                               maxNum = 500, fdrMethod = "BH", sigMethod = "fdr",
                               fdrThr = 0.05, topThr = 10) {
  # INFO: Code is almost identical to oraEnrichment.R, but modified to work with lists.
  #       Additionally, have to get a special intG for the meta-analysis, since it is a merge of the other intG's.
  for (i in seq_along(referenceGene)) {
    referenceGene[[i]] <- intersect(referenceGene[[i]], geneSet[[i]]$gene)
    geneSet[[i]] <- filter(geneSet[[i]], .data$gene %in% referenceGene[[i]])
  }
  geneSetNum <- lapply(geneSet, function(x) {
    tapply(x$gene, x$geneSet, length)
  })

  geneSetNum <- lapply(geneSetNum, function(x) {
    x[x >= minNum & x <= maxNum]
  })

  for (i in seq_along(interestGene)) {
    interestGene[[i]] <- intersect(referenceGene[[i]], intersect(interestGene[[i]], geneSet[[i]]$gene))
    if (length(interestGene[[i]]) == 0) {
      stop(paste0("ERROR: No genes in the interesting list at index ", i, " can annote to any functional category."))
    }
  }

  refG <- list()
  intG <- list()
  met_intG <- data.frame()
  for (i in seq_along(geneSetNum)) {
    refG[[i]] <- data.frame(geneSet = names(geneSetNum[[i]]), size = as.numeric(geneSetNum[[i]]), stringsAsFactors = FALSE)
    intG[[i]] <- filter(geneSet[[i]], .data$gene %in% interestGene[[i]])
    if (i == 1) {
      met_intG <- intG[[i]]
    } else {
      met_intG <- rbind(met_intG, intG[[i]])
    }
  }
  met_intG <- distinct(met_intG)
  met_intG <- tapply(met_intG$gene, met_intG$geneSet, paste, collapse = ";")
  met_intG <- data.frame(geneSet = names(met_intG), overlapId = as.character(met_intG), stringsAsFactors = FALSE)
  intGId <- lapply(intG, function(x) {
    tapply(x$gene, x$geneSet, paste, collapse = ";")
  })
  intGId <- lapply(intGId, function(x) {
    data.frame(geneSet = names(x), overlapId = as.character(x), stringsAsFactors = FALSE)
  })
  geneSetFilter <- lapply(geneSet, function(x) {
    x %>%
      filter(!is.na(.data$geneSet)) %>%
      filter(.data$geneSet %in% names(geneSetNum[[i]])) %>%
      select(.data$geneSet, link = .data$description) %>%
      distinct()
  })
  geneSet <- lapply(geneSet, function(x) {
    x[x$geneSet %in% geneSetFilter[[i]]$geneSet, ]
  })
  genes <- lapply(geneSet, function(x) {
    tapply(x$gene, x$geneSet, rbind)
  })
  modified_geneset <- lapply(genes, function(x) {
    names(x)
  })
  rust_result <- rust_multiomics_ora(modified_geneset, genes, interestGene, referenceGene, "fisher")
  rust_result_df <- lapply(rust_result, function(x) {
    data.frame(
      FDR = p.adjust(x$p, method = fdrMethod), pValue = x$p, expect = x$expect,
      enrichmentRatio = x$enrichment_ratio, geneSet = x$gene_set, overlap = x$overlap
    )
  })
  enrichedResultList <- list()
  backgroundList <- list()
  for (i in seq_along(rust_result_df)) {
    if (i == 1) { # Meta-analysis
      enrichedResult <- rust_result_df[[i]] %>%
        left_join(met_intG, by = "geneSet") %>% # get overlapping gene IDs
        arrange(.data$FDR, .data$pValue, .data$enrichmentRatio)
      enrichedResult$overlap <- sapply(enrichedResult$overlapId, function(x) {
        length(unlist(strsplit(x, ";")))
      })
      if (sigMethod == "fdr") {
        enrichedResultSig <- filter(enrichedResult, .data$FDR < fdrThr)
        if (nrow(enrichedResultSig) == 0) {
          warning("No significant gene set is identified based on FDR ", fdrThr, "!")
          enrichedResultList[[i]] <- NULL
          backgroundList[[i]] <- NULL
        } else {
          enrichedResultInsig <- enrichedResult %>%
            filter(.data$FDR >= fdrThr, .data$overlap != 0) %>%
            select(.data$geneSet, .data$enrichmentRatio, .data$FDR, .data$overlap)
          enrichedResultList[[i]] <- enrichedResultSig
          backgroundList[[i]] <- enrichedResultInsig
        }
      } else {
        # for the top method, we only select the terms with at least one annotated interesting gene
        enrichedResult <- enrichedResult %>% filter(.data$overlap != 0)
        if (nrow(enrichedResult) > topThr) {
          enrichedResultSig <- enrichedResult[1:topThr, ]
          enrichedResultInsig <- enrichedResult[(topThr + 1):nrow(enrichedResult), c("geneSet", "enrichmentRatio", "FDR", "overlap")]
        } else {
          enrichedResultSig <- enrichedResult
          enrichedResultInsig <- data.frame()
        }
        enrichedResultList[[i]] <- enrichedResultSig
        backgroundList[[i]] <- enrichedResultInsig
      }
    } else {
      enrichedResult <- geneSetFilter[[i - 1]] %>%
        left_join(refG[[i - 1]], by = "geneSet") %>%
        left_join(
          rust_result_df[[i - 1]],
          by = "geneSet",
        ) %>%
        left_join(intGId[[i - 1]], by = "geneSet") %>% # get overlapping gene IDs
        arrange(.data$FDR, .data$pValue, .data$enrichmentRatio)
      if (sigMethod == "fdr") {
        enrichedResultSig <- filter(enrichedResult, .data$FDR < fdrThr)
        if (nrow(enrichedResultSig) == 0) {
          warning("No significant gene set is identified based on FDR ", fdrThr, "!")
          enrichedResultList[[i]] <- NULL
          backgroundList[[i]] <- NULL
        } else {
          enrichedResultInsig <- enrichedResult %>%
            filter(.data$FDR >= fdrThr, .data$overlap != 0) %>%
            select(.data$geneSet, .data$enrichmentRatio, .data$FDR, .data$overlap)
          enrichedResultList[[i]] <- enrichedResultSig
          backgroundList[[i]] <- enrichedResultInsig
        }
      } else {
        # for the top method, we only select the terms with at least one annotated interesting gene
        enrichedResult <- enrichedResult %>% filter(.data$overlap != 0)
        if (nrow(enrichedResult) > topThr) {
          enrichedResultSig <- enrichedResult[1:topThr, ]
          enrichedResultInsig <- enrichedResult[(topThr + 1):nrow(enrichedResult), c("geneSet", "enrichmentRatio", "FDR", "overlap")]
        } else {
          enrichedResultSig <- enrichedResult
          enrichedResultInsig <- data.frame()
        }
        enrichedResultList[[i]] <- enrichedResultSig
        backgroundList[[i]] <- enrichedResultInsig
      }
    }
  }
  return(list(enriched = enrichedResultList, background = backgroundList))
}
