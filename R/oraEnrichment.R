#' @importFrom dplyr filter select left_join mutate arrange %>% group_by inner_join
#' @importFrom stats p.adjust phyper
oraEnrichment <- function(interestGene, referenceGene, geneSet, minNum = 10, maxNum = 500, fdrMethod = "BH", sigMethod = "fdr", fdrThr = 0.05, topThr = 10) {
  # before running this code, the main code has checked the overlap among interestGene, referenceGene and geneSet.
  # And this three sets should have overlapping genes.

  # The calculation is based on the genes (input and reference)
  # with annotations in GMT (i.e. effective genes shown in final HTML report)
  # One common question is why input GMT affects results
  # While GSEA does not have this
  referenceGene <- intersect(referenceGene, geneSet$gene)

  geneSet <- filter(geneSet, .data$gene %in% referenceGene)

  geneSetNum <- tapply(geneSet$gene, geneSet$geneSet, length)
  geneSetNum <- geneSetNum[geneSetNum >= minNum & geneSetNum <= maxNum]
  if (length(geneSetNum) == 0) {
    stop("ERROR: The number of annotated genes for all functional categories are not from ", minNum, " to ", maxNum, " for the ORA enrichment method.")
  }

  interestGene <- intersect(interestGene, geneSet$gene)
  interestGene <- intersect(interestGene, referenceGene)
  if (length(interestGene) == 0) {
    stop("ERROR: No genes in the interesting list can annotate to any functional category.")
  }
  refG <- data.frame(geneSet = names(geneSetNum), size = as.numeric(geneSetNum), stringsAsFactors = FALSE)
  intG <- filter(geneSet, .data$gene %in% interestGene)
  intGId <- tapply(intG$gene, intG$geneSet, paste, collapse = ";")
  intGId <- data.frame(geneSet = names(intGId), overlapId = as.character(intGId), stringsAsFactors = FALSE)
  geneSetFilter <- geneSet %>%
    filter(!is.na(.data$geneSet)) %>%
    filter(.data$geneSet %in% names(geneSetNum)) %>%
    select(.data$geneSet, link = .data$description) %>%
    distinct()
  geneSet <- geneSet[geneSet$geneSet %in% geneSetFilter$geneSet, ]
  genes <- tapply(geneSet$gene, geneSet$geneSet, rbind)
  rust_result <- ora_rust(names(genes), genes, interestGene, referenceGene)
  rust_result_df <- data.frame(
    FDR = p.adjust(rust_result$p, method = fdrMethod), pValue = rust_result$p, expect = rust_result$expect,
    enrichmentRatio = rust_result$enrichment_ratio, geneSet = rust_result$gene_set, overlap = rust_result$overlap
  )

  enrichedResult <- geneSetFilter %>%
    left_join(refG, by = "geneSet") %>%
    left_join(
      rust_result_df,
      by = "geneSet",
    ) %>%
    left_join(intGId, by = "geneSet") %>% # get overlapping gene IDs
    arrange(.data$FDR, .data$pValue, .data$enrichmentRatio)
  if (sigMethod == "fdr") {
    enrichedResultSig <- filter(enrichedResult, .data$FDR < fdrThr)
    if (nrow(enrichedResultSig) == 0) {
      warning("No significant gene set is identified based on FDR ", fdrThr, "!")
      return(NULL)
    } else {
      enrichedResultInsig <- enrichedResult %>%
        filter(.data$FDR >= fdrThr, .data$overlap != 0) %>%
        select(.data$geneSet, .data$enrichmentRatio, .data$FDR, .data$overlap)
      return(list(enriched = enrichedResultSig, background = enrichedResultInsig))
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
    return(list(enriched = enrichedResultSig, background = enrichedResultInsig))
  }
}
