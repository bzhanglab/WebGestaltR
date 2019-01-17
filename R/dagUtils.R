#' Expand enriched GO IDs to include ancestors up to the root
#'
#' Returns expanded nodes and DAG tree edges
#'
#' @importFrom dplyr filter
#' @keywords internal
#'
expandDag <- function(goTermList, dagEdgeList) {
	colnames(dagEdgeList) <- c("source", "target")
	queue <- as.list(goTermList)
	## expand to include all linked nodes
	edges <- list()
	while (length(queue) > 0 ) {
		goTerm <- queue[[1]]
		if (length(queue) == 1) {
			queue <- list()
		} else {
			queue <- queue[2:length(queue)]
		}
		inEdges <- filter(dagEdgeList, .data$target == goTerm)
		if (nrow(inEdges) > 0) {
			edges <- c(edges, lapply(split(inEdges, seq(nrow(inEdges))), function(x) list("data"=x)))
		}
		for (parentNode in inEdges$source) {
			if (!parentNode %in% goTermList) {
				queue <- c(queue, parentNode)
				goTermList <- c(goTermList, parentNode)
			}
		}
	}

	return(list(allNodes=goTermList, edges=edges))
}

#' @importFrom dplyr filter
getDagNodes <- function(enrichedRes, allGoList, goIdName, enrichMethod, dagColorSchema) {
	if (!is.null(goIdName)) {
		colnames(goIdName) <- c("id", "name")
	}
	palette <- getColorPalette(enrichedRes, enrichMethod, dagColorSchema)
	return(lapply(allGoList, function(x) {
		goName <- ifelse(is.null(goIdName), "", filter(goIdName, .data$id == x)[[1, "name"]])
		color <- palette(x)
		return(list(
			data=list(
				id=x,
				name=goName,
				color=color
			)
		))
	}))
}


colorPos <- "steelblue"
colorNeg <- "darkorange"
colorNeutral <- "white"

#' @importFrom dplyr filter
getColorPalette <- function(enrichedRes, enrichMethod, schema) {
	if (schema == "binary") {
		if (enrichMethod == "ORA") {
			return(function(goTerm) {
				return(ifelse(goTerm %in% enrichedRes$geneSet, colorPos, colorNeutral))
			})
		} else if (enrichMethod == "GSEA") {
			return(function(goTerm) {
				nes <- filter(enrichedRes, .data$geneSet == goTerm)[["NES"]]
				if (length(nes) == 0) {
					return(colorNeutral)
				} else {
					return(ifelse(nes > 0, colorPos, colorNeg))
				}
			})
		}
	} else if (schema == "continuous") {
		if (enrichMethod=="ORA") {
			minFdr <- min(enrichedRes$FDR)
			minFdrLog <- ifelse(minFdr==0, -log10(.Machine$double.eps), -log10(minFdr))
			colorPalette <- colorRampPalette(c(colorNeutral, colorPos))(128)
			myBreak <- seq(0, minFdrLog + 0.01, length.out=129)

			return(function(goTerm) {
				fdr <- filter(enrichedRes, .data$geneSet == goTerm)[["FDR"]]
				if (length(fdr) == 0) {
					return(colorNeutral)
				} else {
					fdrLog <- ifelse(fdr == 0, -log10(.Machine$double.eps), -log10(fdr))
					return(colorPalette[max(which(myBreak <= fdrLog))])
				}
			})
		} else if (enrichMethod=="GSEA") {
			fdr <- enrichedRes$FDR
			fdr[fdr == 0] <- .Machine$double.eps
			fdr <- sign(enrichedRes$NES) * (-log10(fdr))
			minFdrLog <- min(fdr)
			maxFdrLog <- max(fdr)
			tmp <- getPaletteForGsea(maxFdrLog, minFdrLog)
			colorPalette <- tmp[[1]]
			myBreak <- tmp[[2]]
			return(function(goTerm) {
				row <- filter(enrichedRes, .data$geneSet == goTerm)
				if (nrow(row) == 0) {
					return(colorNeutral)
				} else {
					fdr <- row[["FDR"]]
					nes <- row[["NES"]]
					fdrLog <- ifelse(fdr == 0, sign(nes) * (-log10(.Machine$double.eps)), sign(nes) * (-log10(fdr)))
					return(colorPalette[max(which(myBreak <= fdrLog))])
				}
			})
		}
	}
}

getPaletteForGsea <- function(maxScore, minScore) {
	if (minScore > 0) {
		colorPalette <- colorRampPalette(c(colorNeutral, colorPos))(128)
		myBreak <- seq(0, maxScore + 0.01, length.out=129)
	}else{
		if (maxScore < 0) {
			colorPalette <- colorRampPalette(c(colorNeg, colorNeutral))(128)
			myBreak <- seq(minScore- 0.01, 0, length.out=129)
		}else{
			if (abs(minScore) > maxScore) {
				colorPalette <- colorRampPalette(c(colorNeg, colorNeutral, colorPos))(256)
				myBreak <- c(seq(minScore-0.01, -0.01, length.out=128), 0, seq(0.01, -minScore+0.01, length.out=128))
			}else{
				colorPalette <- colorRampPalette(c(colorNeg, colorNeutral, colorPos))(256)
				myBreak <- c(seq(-maxScore-0.01, -0.01, length.out=128), 0, seq(0.01, maxScore+0.01, length.out=128))
			}
		}
	}
	return(list(colorPalette, myBreak))
}
