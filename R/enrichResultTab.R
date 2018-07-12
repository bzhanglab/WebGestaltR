enrichResultTab <- function(enrichedSig, background) {
	data <- list()
	template <- readLines(system.file("inst/templates/enrichResultTabPlot.mustache", package="WebGestaltR"))
	return(whisker.render(template, data))
}
