mapUserId <- function(enrichedSig,geneColumn,interestingGeneMap){
	####map entrez gene back to the original user id and add one more column to the enrichedSig
	standardId <- interestingGeneMap$standardId
	mapgene <- interestingGeneMap$mapped[,c("userid",standardId)]
	gene <- enrichedSig[,geneColumn]
	gene <- strsplit(gene,";")
	gene <- unlist(lapply(gene,geneM,mapgene))
	enrichedSig <- data.frame(enrichedSig, UserID=gene, stringsAsFactors=FALSE)
	return(enrichedSig)
}

geneM <- function(geneList,mappingTable){
	if(length(geneList)==1 && is.na(geneList)){
		###The categories outputted from GSEA may not have leading edge genes
		return(NA)
	}else{
		u <- mappingTable[mappingTable[,2] %in% geneList,1]
		u <- paste(u,collapse=";")
		return(u)
	}
}

getGeneTables <- function(organism, enrichedSig, geneColumn, interestingGeneMap) {
	if (organism != "others") {
		standardId <- interestingGeneMap$standardId
		mapping <- select(interestingGeneMap$mapped, userid, genesymbol, genename, glink, standardId)
	}
	table <- list()
	for (i in 1:nrow(enrichedSig)) {
		genes <- enrichedSig[[i, geneColumn]]
		geneSetId <- enrichedSig[[i, "geneset"]]
		if (length(genes) == 1 && is.na(genes)) {
			table[[geneSetId]] <- list()
		} else {
			genes <- unlist(strsplit(genes, ";"))
			if (organism != "others") {
				table[[geneSetId]] <- unname(rowSplit(mapping[mapping[[standardId]] %in% genes, ]))
			} else {
				table[[geneSetId]] <- unname(rowSplit(data.frame("userid"=genes)))
			}
		}
	}
	return(table)
}

