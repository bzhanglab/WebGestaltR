#' @importFrom whisker whisker.render
goSlimReport <- function(projectName) {
	goSlimPicPath <- paste0("goslim_summary_", projectName, ".png")
	template <- readLines(system.file("templates/goSlimReport.mustache", package="WebGestaltR"))
	return(whisker.render(template, list(goSlimPicPath=goSlimPicPath)))
}
