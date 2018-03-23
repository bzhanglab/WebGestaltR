detailedGeneTable <- function(outputHtmlFile,organism,enrichMethod,fdrMethod,enrichedSigSub,geneSetDes,geneSetNet,interestingGeneMap,enrichDatabase){

	standardId <- interestingGeneMap$standardId
	nRow <- nrow(enrichedSigSub)
	enrichedCategories <- cbind(enrichedSigSub, "index"=c(1:nRow))

	if(enrichMethod=="ORA") {
		idName <- "overlapID"
		data <- list("C"=enrichedSigSub[,"C"], "O"=enrichedSigSub[,"O"], "E"=round(enrichedSigSub[,"E"],digits=2),
					"R"=round(enrichedSigSub[,"R"],digits=2),
					"PValue"=format(enrichedSigSub[,"PValue"],scientific=TRUE,digits=3),
					"FDR"=format(enrichedSigSub[,"FDR"],scientific=TRUE,digits=3),
					"methodIsOra"=TRUE, "methodIsGsea"=FALSE
					)
	} else if (enrichMethod=="GSEA") {
		idName <- "leadingEdgeID"
		data <- list("Size"=enrichedSigSub[,"Size"], "L"=enrichedSigSub[,"leadingEdgeNum"],
					"ES"=round(enrichedSigSub[,"ES"],digits=2), "NES"=round(enrichedSigSub[,"NES"],digits=2),
					"Pvalue"=format(enrichedSigSub[,"PValue"],scientific=TRUE,digits=3),
					"FDR"=format(enrichedSigSub[,"FDR"],scientific=TRUE,digits=3),
					"methodIsOra"=FALSE, "methodIsGsea"=TRUE
					)
	}
	data$organism <- organism
	data$hasGeneSetDes <- !is.null(geneSetDes)

	data$colNames <- colnames(interestingGeneMap$mapped)[1:4]
	data$numCol <- length(data$colNames)

	enrichedCategories$hasEnrichedGenes <- unname(!is.na(enrichedSigSub[, idName]))
	tableContent <- vector("list", nRow)
	eLinks <- vector("character", nRow)
	for(i in c(1:nRow)){
		g <- enrichedSigSub[i,idName]
		if(!is.na(g)){  ##some enriched terms from GSEA may not have leading edge genes
			g <- unlist(strsplit(g,";"))
			x <- interestingGeneMap$mapped
			names <- colnames(x)
			names[which(names==standardId)] <- 'idCol'
			colnames(x) <- names
			tableContent[[i]] <- unname(rowSplit(x[x[, 'idCol'] %in% g,,drop=FALSE]))
		}
		eLinks[[i]] <- linkModification(enrichDatabase,enrichedSigSub[i,"link"],enrichedSigSub[i,idName],interestingGeneMap)
	}
	enrichedCategories$eLink <- eLinks

	data$hasGeneSetNet <- !is.null(geneSetNet)
	if(data$hasGeneSetNet) {
		enrichedCategories$geneSymbolList <- sapply(tableContent,function(x) {y <- unlist(x); return(paste(y[names(y)=='genesymbol'],collapse=";"))})
		if(enrichMethod=="GSEA"){
			enrichedCategories$scoreList <- sapply(tableContent,function(x) {y <- unlist(x); return(paste(y[names(y)=="score"],collapse=";"))})
		}
	}

	enrichedCategories <- unname(rowSplit(enrichedCategories))
	data$enrichedCategories <- mapply(function(x, y) append(x, list(tableContent=y)), enrichedCategories, tableContent, SIMPLIFY=FALSE)

	statDes <- readLines(system.file("inst/templates/enrichResultStat.mustache", package="WebGestaltR"))
	template <- readLines(system.file("inst/templates/detailedGeneTable.mustache", package="WebGestaltR"))
	return(whisker.render(template, data, partials=list(statDes=statDes)))
}
