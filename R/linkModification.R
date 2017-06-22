linkModification <- function(enrichDatabase,enrichpathwaylink,genelist,interestingGeneMap){
	#####Modify the link to highlight the genes in the pathways. Currently, we only have wikipathway and kegg pathways that need to modify the link########

	if(enrichDatabase=="pathway_KEGG"){
		link <- keggLinkModification(enrichpathwaylink,genelist)
		return(link)
	}
	if(enrichDatabase=="pathway_Wikipathway"){
		link <- wikiLinkModification(enrichpathwaylink,genelist,interestingGeneMap)
		return(link)
	}
	return(enrichpathwaylink)
	
}

keggLinkModification <- function(enrichpathwaylink,genelist){
	genelist <- gsub(";","+",genelist)
	enrichpathwaylink <- paste(enrichpathwaylink,"+",genelist,sep="")
	return(enrichpathwaylink)
}

wikiLinkModification <- function(enrichpathwaylink,genelist,interestingGeneMap){
	genemap <- interestingGeneMap$mapped
	genelist <- unlist(strsplit(genelist,";"))
	gene_symbol <- genemap[genemap[,"entrezgene"] %in% genelist,"genesymbol"]
	for(i in c(1:length(gene_symbol))){
		enrichpathwaylink <- paste(enrichpathwaylink,"&label[]=",gene_symbol[i],sep="")
	}
	enrichpathwaylink <- paste(enrichpathwaylink,"&colors=red",sep="")
	return(enrichpathwaylink)
}
