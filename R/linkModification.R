linkModification <- function(enrichDatabase,enrichPathwayLink,geneList,interestingGeneMap){
	#####Modify the link to highlight the genes in the pathways. Currently, we only have wikipathway and kegg pathways that need to modify the link########

	if(enrichDatabase=="pathway_KEGG"){
		link <- keggLinkModification(enrichPathwayLink,geneList)
		return(link)
	}
	if(startsWith(enrichDatabase, "pathway_Wikipathway")){
		link <- wikiLinkModification(enrichPathwayLink,geneList,interestingGeneMap)
		return(link)
	}
	return(enrichPathwayLink)
}

keggLinkModification <- function(enrichPathwayLink,geneList){
	geneList <- gsub(";","+",geneList)
	enrichPathwayLink <- paste(enrichPathwayLink,"+",geneList,sep="")
	return(enrichPathwayLink)
}

wikiLinkModification <- function(enrichPathwayLink,geneList,interestingGeneMap){
	geneMap <- interestingGeneMap$mapped
	geneList <- unlist(strsplit(geneList,";"))
	for(i in c(1:length(geneList))) {
		enrichPathwayLink <- paste0(enrichPathwayLink, "&xref[]=", geneList[i], ",Entrez Gene")
	}
	enrichPathwayLink <- paste0(enrichPathwayLink, "&colors=red")
	return(enrichPathwayLink)
}
