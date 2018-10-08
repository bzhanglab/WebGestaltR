enrichResultSection <- function(enrichMethod, geneSetDes, geneSetDag, clusters) {
	data <- list(methodIsOra=enrichMethod=='ORA',
				hasGeneSetDes=!is.null(geneSetDes),
				hasGeneSetDag=!is.null(geneSetDag),
				hasAp=!is.null(clusters$ap),
				hasWsc=!is.null(clusters$wsc)
				)
	template <- readLines(system.file("templates/enrichResultSection.mustache", package="WebGestaltR"))
	return(whisker.render(template, data))
}
