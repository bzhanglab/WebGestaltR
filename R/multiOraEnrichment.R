multiOraEnrichment <- function(interestGene, referenceGene, geneSet, minNum = 10, maxNum = 500, fdrMethod = "BH", sigMethod = "fdr", fdrThr = 0.05, topThr = 10) {
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
  }
}
