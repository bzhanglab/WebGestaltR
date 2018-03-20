detailedGeneTable <- function(outputHtmlFile,organism,enrichMethod,fdrMethod,enrichedSig_sub,geneSetDes,geneSetNet,interestingGeneMap,enrichDatabase){

	standardId <- interestingGeneMap$standardId
	nRow <- nrow(enrichedSig_sub)
	enrichedCategories <- cbind(enrichedSig_sub, "index"=c(1:nRow))

	if(enrichMethod=="ORA") {
		idName <- "overlapID"
		data <- list("C"=enrichedSig_sub[,"C"], "O"=enrichedSig_sub[,"O"], "E"=round(enrichedSig_sub[,"E"],digits=2),
					"R"=round(enrichedSig_sub[,"R"],digits=2),
					"PValue"=format(enrichedSig_sub[,"PValue"],scientific=TRUE,digits=3),
					"FDR"=format(enrichedSig_sub[,"FDR"],scientific=TRUE,digits=3),
					"methodIsOra"=TRUE, "methodIsGsea"=FALSE
					)
	} else if (enrichMethod=="GSEA") {
		idName <- "leadingEdgeID"
		data <- list("Size"=enrichedSig_sub[,"Size"], "L"=enrichedSig_sub[,"leadingEdgeNum"],
					"ES"=round(enrichedSig_sub[,"ES"],digits=2), "NES"=round(enrichedSig_sub[,"NES"],digits=2),
					"Pvalue"=format(enrichedSig_sub[,"PValue"],scientific=TRUE,digits=3),
					"FDR"=format(enrichedSig_sub[,"FDR"],scientific=TRUE,digits=3),
					"methodIsOra"=FALSE, "methodIsGsea"=TRUE
					)
	}
	data$organism <- organism
	data$hasGeneSetDes <- !is.null(geneSetDes)

	data$colNames <- colnames(interestingGeneMap$mapped)[1:4]
	data$numCol <- length(data$colNames)

	enrichedCategories$hasEnrichedGenes <- unname(!is.na(enrichedSig_sub[, idName]))
	tableContent <- vector("list", nRow)
	eLinks <- vector("character", nRow)
	for(i in c(1:nRow)){
		g <- enrichedSig_sub[i,idName]
		if(!is.na(g)){  ##some enriched terms from GSEA may not have leading edge genes
			g <- unlist(strsplit(g,";"))
			x <- interestingGeneMap$mapped
			names <- colnames(x)
			names[which(names==standardId)] <- 'idCol'
			colnames(x) <- names
			tableContent[[i]] <- unname(rowSplit(x[x[, 'idCol'] %in% g,,drop=FALSE]))
		}
		eLinks[[i]] <- linkModification(enrichDatabase,enrichedSig_sub[i,"link"],enrichedSig_sub[i,idName],interestingGeneMap)
	}
	enrichedCategories$eLink <- eLinks

	data$hasGeneSetNet <- !is.null(geneSetNet)
	if(data$hasGeneSetNet) {
		enrichedCategories$geneSymbolList <- sapply(tableContent,function(x) {return(paste(x[,"genesymbol"],collapse=";"))})
		if(enrichMethod=="GSEA"){
			enrichedCategories$scoreList <- sapply(tableContent,function(x) {return(paste(x[,"score"],collapse=";"))})
		}
	}

	enrichedCategories <- unname(rowSplit(enrichedCategories))
	data$enrichedCategories <- mapply(function(x, y) append(x, list(tableContent=y)), enrichedCategories, tableContent, SIMPLIFY=FALSE)

	statDes <- readLines(system.file("inst/templates/enrichResultStat.mustache", package="WebGestaltR"))
	template <- readLines(system.file("inst/templates/detailedGeneTable.mustache", package="WebGestaltR"))
	return(whisker.render(template, data, partials=list(statDes=statDes)))
}
