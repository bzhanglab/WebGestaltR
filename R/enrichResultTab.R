enrichResultTab <- function(enrichMethod, geneSetDes, geneSetDag) {
	data <- list(methodIsOra=enrichMethod=='ORA',
				hasGeneSetDes=!is.null(geneSetDes),
				hasGeneSetDag=!is.null(geneSetDag)
				)
	template <- readLines(system.file("templates/enrichResultTabPlot.mustache", package="WebGestaltR"))
	return(whisker.render(template, data))
}
