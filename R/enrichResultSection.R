#' enrichResultSection
#'
#' Conditionally render template of main result section. Actual work is carried out in front end
#'
#' @importFrom whisker whisker.render
#' @keywords internal
#'
enrichResultSection <- function(enrichMethod, geneSetDes, geneSetDag, geneSetNet, clusters) {
	data <- list(methodIsOra=enrichMethod=='ORA',
				hasGeneSetDes=!is.null(geneSetDes),
				hasGeneSetDag=!is.null(geneSetDag),
				hasGeneSetNet=!is.null(geneSetNet),
				hasAp=!is.null(clusters$ap),
				hasWsc=!is.null(clusters$wsc)
				)
	template <- readLines(system.file("templates/enrichResultSection.mustache", package="WebGestaltR"))
	return(whisker.render(template, data))
}
