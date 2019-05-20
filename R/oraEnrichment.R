#' @importFrom dplyr filter select left_join mutate arrange %>%
#' @importFrom stats p.adjust phyper
oraEnrichment <- function(interestGene, referenceGene, geneSet, minNum=10, maxNum=500, fdrMethod="BH", sigMethod="fdr", fdrThr=0.05, topThr=10) {
	#before running this code, the main code has checked the overlap among interestGene, referenceGene and geneSet.
	#And this three sets should have overlapping genes.

	referenceGene <- intersect(referenceGene, geneSet$gene)

	geneSet <- filter(geneSet, .data$gene %in% referenceGene)

	geneSetNum <- tapply(geneSet$gene, geneSet$geneSet,length)
	geneSetNum <- geneSetNum[geneSetNum>=minNum & geneSetNum<=maxNum]
	if (length(geneSetNum) == 0) {
		stop("ERROR: The number of annotated genes for all functional categories are not from ", minNum, " to ", maxNum, " for the ORA enrichment method.")
	}

	geneSet <-filter(geneSet, .data$geneSet %in% names(geneSetNum))

	interestGene <- intersect(interestGene, geneSet$gene)
	interestGene <- intersect(interestGene, referenceGene)
	if (length(interestGene) == 0) {
		stop("ERROR: No genes in the interesting list can annotate to any functional category.")
	}

	###############Enrichment analysis###################
	ra <- length(interestGene) / length(referenceGene)
	refG <- data.frame(geneSet=names(geneSetNum), size=unname(geneSetNum), stringsAsFactors=FALSE)

	intG <- filter(geneSet, .data$gene %in% interestGene)
	intGNum <- tapply(intG$gene, intG$geneSet, length)
	intGNum <- data.frame(geneSet=names(intGNum), overlap=unname(intGNum), stringsAsFactors=FALSE)

	intGId <- tapply(intG$gene, intG$geneSet, paste, collapse=";")
	intGId <- data.frame(geneSet=names(intGId), overlapId=unname(intGId), stringsAsFactors=FALSE)

	enrichedResult <- geneSet %>% filter(!is.na(.data$geneSet)) %>%
		select(.data$geneSet, link=.data$description) %>% distinct() %>%
		left_join(refG, by="geneSet") %>%
		left_join(intGNum, by="geneSet") %>% # this may just be inner_join. O is NA should not be meaningful anyway
		mutate(overlap=ifelse(is.na(.data$overlap), 0, .data$overlap), expect=.data$size * ra, enrichmentRatio=.data$overlap /.data$expect,
			   pValue=1-phyper(.data$overlap - 1, length(interestGene), length(referenceGene) - length(interestGene), .data$size, lower.tail=TRUE, log.p=FALSE),
			   FDR=p.adjust(.data$pValue, method=fdrMethod)
		) %>%
		left_join(intGId, by="geneSet") %>%
		arrange(.data$FDR, .data$pValue)

	if (sigMethod == "fdr") {
		enrichedResultSig <- filter(enrichedResult, .data$FDR<fdrThr)
		if (nrow(enrichedResultSig) == 0) {
			warning("No significant gene set is identified based on FDR ", fdrThr, "!")
			return(NULL)
		} else {
			enrichedResultInsig <- enrichedResult %>% filter(.data$FDR >= fdrThr, .data$overlap != 0) %>% select(.data$geneSet, .data$enrichmentRatio, .data$FDR, .data$overlap)
			return(list(enriched=enrichedResultSig, background=enrichedResultInsig))
		}
	} else {
		#for the top method, we only select the terms with at least one annotated interesting gene
		enrichedResult <- enrichedResult %>% filter(.data$overlap != 0)
		if (nrow(enrichedResult)>topThr) {
			enrichedResultSig <- enrichedResult[1:topThr, ]
			enrichedResultInsig <- enrichedResult[(topThr+1):nrow(enrichedResult), c("geneSet", "enrichmentRatio", "FDR", "overlap")]
		}else{
			enrichedResultSig <- enrichedResult
			enrichedResultInsig <- data.frame()
		}
		return(list(enriched=enrichedResultSig, background=enrichedResultInsig))
	}
}
