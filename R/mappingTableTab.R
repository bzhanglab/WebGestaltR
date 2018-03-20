mappingTableTab <- function(interestingGeneMap){
	standardId <- interestingGeneMap$standardId
	tableNames <- colnames(interestingGeneMap$mapped)[1:4]
	colnames(interestingGeneMap$mapped)[4] <- 'idCol'  # temporary change for template rendering
	mappedGenes <- unname(rowSplit(interestingGeneMap$mapped))
	unmappedGenes <- interestingGeneMap$unmapped

	template <- readLines(system.file("inst/templates/mappingTableTab.mustache", package="WebGestaltR"))
	data <- list(tableNames=tableNames, mappedGenes=mappedGenes, unmappedGenes=unmappedGenes, standardId=standardId)
	table <- whisker.render(template, data=data)
	colnames(interestingGeneMap$mapped)[4] <- standardId

	return(table)
}
