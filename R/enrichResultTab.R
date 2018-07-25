enrichResultTab <- function(enrichMethod) {
	data <- list(methodIsOra=enrichMethod=='ORA')
	template <- readLines(system.file("inst/templates/enrichResultTabPlot.mustache", package="WebGestaltR"))
	return(whisker.render(template, data))
}
