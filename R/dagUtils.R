expandDag <- function(goTermList, dagEdgeList) {
## Expand enriched GO IDs to include ancestors up to the root
## returns expanded nodes and DAG tree edges.
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
		inEdges <- filter(dagEdgeList, target == goTerm)
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

getDagNodes <- function(enrichedRes, allGoList, goIdName, enrichMethod, dagColorSchema) {
	if (!is.null(goIdName)) {
		colnames(goIdName) <- c("id", "name")
	}
	return(lapply(allGoList, function(x) {
		goName <- ifelse(is.null(goIdName), "", filter(goIdName, id == x)[[1, "name"]])
		color <- getDagNodeColor(enrichedRes, x, enrichMethod, dagColorSchema)
		return(list(
			data=list(
				id=x,
				name=goName,
				color=color
			)
		))
	}))
}


getDagNodeColor <- function(enrichedRes, goTerm, enrichMethod, schema) {
	if (schema == "binary") {
		if (enrichMethod == "ORA") {
			return(ifelse(goTerm %in% enrichedRes$geneSet, "red", "white"))
		} else if (enrichMethod == "GSEA") {
			nes <- filter(enrichedRes, geneSet == goTerm)[["NES"]]
			return(ifelse(nes > 0, "red", "blue"))
		}
	} else if (schema == "continuous") {
		if (enrichMethod=="ORA") {
			minFdr <- min(enrichedRes$FDR)
			minFdrLog <- ifelse(minFdr==0, -log10(2.2e-16), -log10(minFdr))
			colorPalette <- colorRampPalette(c("white", "red"))(128)
			myBreak <- seq(0, minFdrLog + 0.01, length.out=129)

			fdr <- filter(enrichedRes, geneSet == goTerm)[["FDR"]]
			if (length(fdr) == 0) {
				return(colorPalette[1])
			}
			fdrLog <- ifelse(fdr == 0, -log10(2.2e-16), -log10(fdr))
			return(colorPalette[max(which(myBreak <= fdrLog))])
		} else if (enrichMethod=="GSEA") {
			fdr <- enrichedRes$FDR
			fdr[fdr == 0] <- 2.2e-16
			fdr <- sign(enrichedRes$NES) * (-log10(fdr))
			minFdrLog <- min(fdr)
			maxFdrLog <- max(fdr)
			if (minFdrLog > 0) {
				colorPalette <- colorRampPalette(c("white", "red"))(128)
				myBreak <- seq(0, maxFdrLog + 0.01, length.out=129)
			}else{
				if (maxFdrLog < 0) {
					colorPalette <- colorRampPalette(c("blue", "white"))(128)
					myBreak <- seq(minFdrLog- 0.01, 0, length.out=129)
				}else{
					if (abs(minFdrLog) > maxFdrLog) {
						colorPalette <- colorRampPalette(c("blue", "white", "red"))(256)
						myBreak <- c(seq(minFdrLog-0.01, -0.01, length.out=128), 0, seq(0.01, -minFdrLog+0.01, length.out=128))
					}else{
						colorPalette <- colorRampPalette(c("blue", "white", "red"))(256)
						myBreak <- c(seq(-maxFdrLog-0.01, -0.01, length.out=128), 0, seq(0.01, maxFdrLog+0.01, length.out=128))
					}
				}
			}
			fdr <- filter(enrichedRes, geneSet == goTerm)[["FDR"]]
			if (length(fdr) == 0) {
				return("#ffffff")
			}
			nes <- filter(enrichedRes, geneSet == goTerm)[["NES"]]
			fdrLog <- ifelse(fdr == 0, sign(nes) * (-log10(2.2e-16)), sign(nes) * (-log10(fdr)))
			return(colorPalette[max(which(myBreak <= fdrLog))])
		}
	}

}
