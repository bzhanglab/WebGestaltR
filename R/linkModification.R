linkModification <- function(enrichDatabase,enrichPathwayLink,geneList,interestingGeneMap){
	#####Modify the link to highlight the genes in the pathways. Currently, we only have wikipathway and kegg pathways that need to modify the link########

	if(enrichDatabase=="pathway_KEGG"){
		link <- keggLinkModification(enrichPathwayLink,geneList)
		return(link)
	}
	if(enrichDatabase=="pathway_Wikipathway"){
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
	geneSymbol <- geneMap[geneMap[,"entrezgene"] %in% geneList,"genesymbol"]
	for(i in c(1:length(geneSymbol))){
		enrichPathwayLink <- paste(enrichPathwayLink,"&label[]=",geneSymbol[i],sep="")
	}
	enrichPathwayLink <- paste(enrichPathwayLink,"&colors=red",sep="")
	return(enrichPathwayLink)
}
