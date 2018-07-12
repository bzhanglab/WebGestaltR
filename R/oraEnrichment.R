oraEnrichment <- function(interestGene,referenceGene,geneSet,minNum=10,maxNum=500,fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,topThr=10){
	#before running this code, the main code has checked the overlap among interestGene, referenceGene and geneSet.
	#And this three sets should have overlapping genes.

	interestGene <- as.character(interestGene)
	referenceGene <- as.character(referenceGene)

	geneSet[,3] <- as.character(geneSet[,3])

	referenceGene <- intersect(referenceGene,geneSet[,3])

	geneSet <- geneSet[geneSet[,3] %in% referenceGene,,drop=FALSE]

	geneSetNum <- tapply(geneSet[,3],geneSet[,1],length)
	geneSetNum <- geneSetNum[geneSetNum>=minNum & geneSetNum<=maxNum]
	if(length(geneSetNum)==0){
		error <- paste("ERROR: The number of annotated genes for all functional categories are not from ",minNum," to ",maxNum," for the ORA enrichment method.",sep="")
		cat(error)
		return(error)
	}

	geneSet <- geneSet[geneSet[,1] %in% names(geneSetNum),,drop=FALSE]

	interestGene <- intersect(interestGene,geneSet[,3])
	interestGene <- intersect(interestGene,referenceGene)

	if(length(interestGene)==0){
		error <- "ERROR: No genes in the interesting list can annotate to any functional category."
		cat(error)
		return(error)
	}


	###############Enrichment analysis###################
	ra <- length(interestGene)/length(referenceGene)

	g <- unique(geneSet[!is.na(geneSet[,1]), c(1,2)])
	enrichedResult <- data.frame(geneset=g[,1],link=g[,2],stringsAsFactors=FALSE)

	refG <- data.frame(geneset=names(geneSetNum),C=unname(geneSetNum),stringsAsFactors=FALSE)
	enrichedResult <- merge(enrichedResult,refG,by="geneset",all.x=TRUE)

	intG <- geneSet[geneSet[,3] %in% interestGene,,drop=FALSE]
	intGNum <- tapply(intG[,3],intG[,1],length)
	intGNum <- data.frame(geneset=names(intGNum),O=unname(intGNum),stringsAsFactors=FALSE)
	enrichedResult <- merge(enrichedResult,intGNum,by="geneset",all.x=TRUE)
	enrichedResult[is.na(enrichedResult[,4]),4] <- 0

	enrichedResult[,5] <- enrichedResult[,3]*ra
	enrichedResult[,6] <- enrichedResult[,4]/enrichedResult[,5]
	enrichedResult[is.na(enrichedResult[,6]),6] <- NA
	enrichedResult[,7] <- 1-phyper(enrichedResult[,4]-1,length(interestGene),length(referenceGene)-length(interestGene),enrichedResult[,3],lower.tail = TRUE,log.p= FALSE)
	enrichedResult[,8] <- p.adjust(enrichedResult[,7],method=fdrMethod)
	colnames(enrichedResult)[5:8] <- c("E","R","PValue","FDR")

	intG <- tapply(intG[,3],intG[,1],paste,collapse=";")
	intG <- data.frame(geneset=names(intG),overlapID=unname(intG),stringsAsFactors=FALSE)
	enrichedResult <- merge(enrichedResult,intG,by="geneset",all.x=TRUE)
	enrichedResult[is.na(enrichedResult[,9]),9] <- NA

	enrichedResult <- enrichedResult[order(enrichedResult[,8],enrichedResult[,7]),]
	if(sigMethod=="fdr"){
		enrichedResultSig <- enrichedResult[enrichedResult[,8]<fdrThr,]
		if(nrow(enrichedResultSig)==0){
			cat("No significant gene set is identified based on FDR ",fdrThr,"!",sep="")
			return(NULL)
		}else{
<<<<<<< HEAD
			enrichedResultInsig <- enrichedResult[enrichedResult[,8]>=fdrThr&enrichedResult$R!=0, c("geneset", "R", "FDR", "O")]
=======
			enrichedResultInsig <- enrichedResult[enrichedResult[,8]>=fdrThr&d$R!=0, c("geneset", "R", "FDR", "O")]
>>>>>>> a8ddbc3030372e4116183949e61e9004c4f5ebda
			return(list(enriched=enrichedResultSig, background=enrichedResultInsig))
		}
	}else{
		#for the top method, we only select the terms with at least one annotated interesting gene
		x <- enrichedResult[!is.na(enrichedResult[,9]),]
		x <- x[order(x[,8],x[,7]),]
		if(nrow(x)>topThr){
			enrichedResultSig <- x[1:topThr,]
			enrichedResultInsig <- x[(topThr+1):nrow(x), c("geneset", "R", "FDR", "O")]
		}else{
			enrichedResultSig <- x
			enrichedResultInsig <- data.frame()
		}
		return(list(enriched=enrichedResultSig, background=enrichedResultInsig))
	}
}
