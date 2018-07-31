enrichResultTab <- function(enrichMethod, geneSetDag) {
	data <- list(methodIsOra=enrichMethod=='ORA',
				hasGeneSetDag=!is.null(geneSetDag)
				)
	template <- readLines(system.file("inst/templates/enrichResultTabPlot.mustache", package="WebGestaltR"))
	return(whisker.render(template, data))
}
