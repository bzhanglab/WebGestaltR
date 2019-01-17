#' Modify the link to highlight the genes in the pathways
#'
#' Currently, we only have wikipathway and kegg pathways that need to modify the link
#'
#' @keywords internal
linkModification <- function(enrichMethod, enrichDatabase, enrichPathwayLink, geneList, interestingGeneMap) {

	if(enrichDatabase=="pathway_KEGG"){
		link <- keggLinkModification(enrichPathwayLink,geneList)
		return(link)
	}
	if(startsWith(enrichDatabase, "pathway_Wikipathway")){
		link <- wikiLinkModification(enrichMethod, enrichPathwayLink, geneList, interestingGeneMap)
		return(link)
	}
	return(enrichPathwayLink)
}

keggLinkModification <- function(enrichPathwayLink,geneList){
	geneList <- gsub(";","+",geneList)
	enrichPathwayLink <- paste(enrichPathwayLink,"+",geneList,sep="")
	return(enrichPathwayLink)
}

wikiLinkModification <- function(enrichMethod, enrichPathwayLink, geneList, interestingGeneMap) {
	geneMap <- interestingGeneMap$mapped
	geneList <- unlist(strsplit(geneList,";"))
	geneMap <- filter(geneMap, .data$entrezgene %in% geneList)
	enrichPathwayLink <- paste0(enrichPathwayLink,
		paste0(sapply(geneMap$geneSymbol, function(x) paste0("&label[]=", x)), collapse="")
		#not many pathway have entrezgene xref. Using both also seem to interfere with coloring
		#paste0(sapply(geneMap$entrezgene, function(x) paste0("&xref[]=", x, ",Entrez Gene")), collapse="")
	)
	if (enrichMethod == "ORA") {
		enrichPathwayLink <- paste0(enrichPathwayLink, "&colors=", colorPos)
	} else if (enrichMethod == "GSEA") {
		scores <- filter(interestingGeneMap$mapped, .data$entrezgene %in% geneList)[["score"]]
		maxScore <- max(scores)
		minScore <- min(scores)
		tmp <- getPaletteForGsea(maxScore, minScore)
		palette <- tmp[[1]]
		breaks <- tmp[[2]]
		colors <- sapply(scores, function(s) palette[max(which(breaks <= s))])
		colorStr <- paste(gsub("#", "%23", colors, fixed=TRUE), collapse=",")
		enrichPathwayLink <- paste0(enrichPathwayLink, "&colors=", colorStr)
	}
	return(enrichPathwayLink)
}
