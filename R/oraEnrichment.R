oraEnrichment <- function(interestGene,referenceGene,geneSet,minNum=10,maxNum=500,fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,topThr=10){
	#before running this code, the main code has checked the overlap among interestGene, referenceGene and geneSet.
	#And this three sets should have overlapping genes.

	referenceGene <- intersect(referenceGene, geneSet$gene)

	geneSet <- filter(geneSet, gene %in% referenceGene)

	geneSetNum <- tapply(geneSet$gene, geneSet$geneSet,length)
	geneSetNum <- geneSetNum[geneSetNum>=minNum & geneSetNum<=maxNum]
	if(length(geneSetNum)==0){
		error <- paste("ERROR: The number of annotated genes for all functional categories are not from ",minNum," to ",maxNum," for the ORA enrichment method.",sep="")
		cat(error)
		return(error)
	}

	geneSet <-filter(geneSet, geneSet %in% names(geneSetNum))

	interestGene <- intersect(interestGene, geneSet$gene)
	interestGene <- intersect(interestGene, referenceGene)

	if(length(interestGene)==0){
		error <- "ERROR: No genes in the interesting list can annotate to any functional category."
		cat(error)
		return(error)
	}


	###############Enrichment analysis###################
	ra <- length(interestGene)/length(referenceGene)
	refG <- data.frame(geneset=names(geneSetNum),C=unname(geneSetNum),stringsAsFactors=FALSE)

	intG <- filter(geneSet, gene %in% interestGene)
	intGNum <- tapply(intG$gene, intG$geneSet, length)
	intGNum <- data.frame(geneset=names(intGNum),O=unname(intGNum),stringsAsFactors=FALSE)

	intGId <- tapply(intG$gene, intG$geneSet, paste, collapse=";")
	intGId <- data.frame(geneset=names(intGId), overlapID=unname(intGId), stringsAsFactors=FALSE)

	enrichedResult <- geneSet %>% filter(!is.na(geneSet)) %>%
		select(geneset=geneSet, link=description) %>% distinct() %>%
		left_join(refG, by="geneset") %>%
		left_join(intGNum, by="geneset") %>% # this may just be inner_join. O is NA should not be meaningful anyway
		mutate(O=ifelse(is.na(O), 0, O), E=C*ra, R=O/E,
			   PValue=1-phyper(O-1, length(interestGene), length(referenceGene)-length(interestGene), C, lower.tail=TRUE, log.p=FALSE),
			   FDR=p.adjust(PValue, method=fdrMethod)
		) %>%
		left_join(intGId, by="geneset") %>%
		arrange(FDR, PValue)

	if(sigMethod=="fdr"){
		enrichedResultSig <- filter(enrichedResult, FDR<fdrThr)
		if(nrow(enrichedResultSig)==0){
			cat("No significant gene set is identified based on FDR ",fdrThr,"!",sep="")
			return(NULL)
		}else{
			enrichedResultInsig <- enrichedResult %>% filter(FDR>=fdrThr, O!=0) %>% select(geneset, R, FDR, O)
			return(list(enriched=enrichedResultSig, background=enrichedResultInsig))
		}
	}else{
		#for the top method, we only select the terms with at least one annotated interesting gene
		enrichedResult <- enrichedResult %>% filter(O!=0)
		if (nrow(enrichedResult)>topThr) {
			enrichedResultSig <- enrichedResult[1:topThr, ]
			enrichedResultInsig <- enrichedResult[(topThr+1):nrow(enrichedResult), c("geneset", "R", "FDR", "O")]
		}else{
			enrichedResultSig <- enrichedResult
			enrichedResultInsig <- data.frame()
		}
		return(list(enriched=enrichedResultSig, background=enrichedResultInsig))
	}
}
