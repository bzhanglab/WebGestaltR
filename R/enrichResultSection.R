#' enrichResultSection
#'
#' Conditionally render template of main result section. Actual work is carried out in front end
#'
#' @importFrom whisker whisker.render
#' @importFrom dplyr select distinct filter
#' @importFrom jsonlite toJSON
#' @keywords internal
#'
enrichResultSection <- function(enrichMethod, enrichedSig, geneSet, geneSetDes, geneSetDag, geneSetNet, clusters) {
	if ('database' %in% colnames(geneSet)) {
		#multiple databases
		netDatabases <- names(geneSetNet[!sapply(geneSetNet, is.null)])
		setSource <- geneSet %>% select(.data$geneSet, .data$database) %>%
			distinct() %>%
			filter(.data$geneSet %in% enrichedSig$geneSet)
		setsWithNetJson <- toJSON((filter(setSource, .data$database %in% netDatabases))$geneSet)
		hasGeneSetDag <- length(geneSetDag[!sapply(geneSetDag, is.null)]) > 0
		hasMultipleDatabases <- TRUE
	} else {
		setsWithNetJson <- toJSON(!is.null(geneSetNet), auto_unbox=TRUE)
		hasGeneSetDag <- !is.null(geneSetDag)
		hasMultipleDatabases <- FALSE
	}

	data <- list(methodIsOra=enrichMethod=='ORA',
				hasGeneSetDag=hasGeneSetDag,
				hasMultipleDatabases=hasMultipleDatabases,
				setsWithNetJson=setsWithNetJson,
				hasAp=!is.null(clusters$ap),
				hasWsc=!is.null(clusters$wsc)
				)
	template <- readLines(system.file("templates/enrichResultSection.mustache", package="WebGestaltR"))
	return(whisker.render(template, data))
}
