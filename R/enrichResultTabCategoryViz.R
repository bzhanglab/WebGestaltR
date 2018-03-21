enrichResultTabCategoryViz <- function(outputHtmlFile,organism,enrichMethod,fdrMethod,enrichedSig,dNum,geneSetDag,geneSetDes,geneSetNet,outputDirectory,timeStamp,dagColor,hostName,interestingGeneMap,enrichDatabase){
	enrichedSigSub <- extractSubSig(enrichMethod,enrichedSig,dNum)  ###extract dNum significant categories. extractSubSig contains hard code

	data <- list(hasGeneSetDag=!is.null(geneSetDag), methodIsGsea=enrichMethod=="GSEA",
				dagColocrIsBinary=dagColor=="binary", standardId=interestingGeneMap$standardId,
				hasGeneSetDes=!is.null(geneSetDes), methodIsOra=enrichMethod=="ORA"
				)

	dagInfo <- NULL
	enrichedCategories <- NULL
	if(!is.null(geneSetDag)){
		dagInfo <- goDagViz(enrichedSigSub,enrichMethod,geneSetDes,geneSetDag,outputHtmlFile,outputDirectory,timeStamp,dagColor,hostName)
		if(!is.null(enrichedSig) && !is.null(geneSetDag)){
			data$jId <- paste("Project_",timeStamp,"_.json",sep="")
			data$jF <- dagInfo$jF
			data$style <- dagInfo$style
		}
		data$enrichedCategories <- enrichedCategories
	}else{
		#########Create table to summary enriched results###########
		if(enrichMethod=="ORA"){
			enrichedCategories <- cbind(enrichedSigSub, "color"=rep("black", nrow(enrichedSigSub)))
		} else if(enrichMethod=="GSEA"){
			color <- vector("character", nrow(enrichedSigSub))
			color[enrichedSigSub[,"NES"]>0] <- "red"
			color[enrichedSigSub[,"NES"]<=0] <- "blue"
			enrichedCategories <- cbind(enrichedSigSub, "color"=color)
		}
		enrichedCategories[,"FDR"] <- format(enrichedCategories[, "FDR"], scientific=TRUE, digits=3)
		data$enrichedCategories <- unname(rowSplit(enrichedCategories))
	}

	data$detailedGeneTableContent <- detailedGeneTable(outputHtmlFile,organism,enrichMethod,fdrMethod,enrichedSigSub,geneSetDes,geneSetNet,interestingGeneMap,enrichDatabase)

	template <- readLines(system.file("inst/templates/enrichResultTab.mustache", package="WebGestaltR"))
	return(whisker.render(template, data))
}

extractSubSig <- function(enrichMethod,enrichedSig,dNum){
	###extract the dNum siganificant categories from the enrichedSig
	###Hard Code######
	if(enrichMethod=="ORA"){
		enrichedSig <- enrichedSig[order(enrichedSig[,"FDR"],enrichedSig[,"PValue"]),]
		if(nrow(enrichedSig)>dNum){
			enrichedSigSub <- enrichedSig[1:dNum,]
		}else{
			enrichedSigSub <- enrichedSig
		}
	}
	if(enrichMethod=="GSEA"){
		x <- enrichedSig[enrichedSig[,"NES"]>0,]
		y <- enrichedSig[enrichedSig[,"NES"]<0,]
		if(nrow(x)>0){
			x <- x[order(x[,"FDR"],x[,"PValue"]),]
			if(nrow(x)>dNum){
				x <- x[1:dNum,]
			}
		}

		if(nrow(y)>0){
			y <- y[order(y[,"FDR"],y[,"PValue"]),]
			if(nrow(y)>dNum){
				y <- y[1:dNum,]
			}
		}

		enrichedSigSub <- rbind(x,y)
	}
	return(enrichedSigSub)
}


enrichResultOthers <- function(enrichMethod,enrichedSig,geneSetDes,fdrMethod,dNum){
	enrichedSigSub <- extractSubSig(enrichMethod,enrichedSig,dNum)
	extractSig <- data.frame(id=enrichedSigSub[,"geneset"],link=enrichedSigSub[,"link"],stringsAsFactors=FALSE)

	if(enrichMethod=="ORA"){
		extractSig$statistic <- paste("C=",enrichedSigSub[,"C"],"; O=",enrichedSigSub[,"O"],"; E=",round(enrichedSigSub[,"E"],digits=2),"; R=",round(enrichedSigSub[,"R"],digits=2),"; PValue=",format(enrichedSigSub[,"P Value"],scientific=TRUE,digits=3),"; FDR=",format(enrichedSigSub[,"FDR"],scientific=TRUE,digits=3),sep="")
		extractSig$genes <- gsub(",|;"," ",enrichedSigSub[,"overlapID"])
		extractSig$color <- rep("black", nrow(extractSig))
	} else if (enrichMethod=="GSEA"){
		extractSig$statistic <- paste("Size=",enrichedSigSub[,"Size"],"; L=",enrichedSigSub[,"leadingEdgeNum"],"; ES=",round(enrichedSigSub[,"ES"],digits=2),"; NES=",round(enrichedSigSub[,"NES"],digits=2),"; P value=",format(enrichedSigSub[,"PValue"],scientific=TRUE,digits=3),"; FDR=",format(enrichedSigSub[,"FDR"],scientific=TRUE,digits=3),sep="")
		extractSig$genes <- gsub(",|;"," ",enrichedSigSub[,"leadingEdgeID"])
		extractSig$color <- sapply(enrichedSigSub[,"NES"], function(x) {if(x>0) {return("red")} else {return("blue")}})
	}

	if(!is.null(geneSetDes)){
		extractSig$name <- enrichedSigSub[,"description"]
	}

	template <- readLines(system.file("inst/templates/enrichResultOthers.mustache", package="WebGestaltR"))
	statDes <- readLines(system.file("inst/templates/enrichResultStat.mustache", package="WebGestaltR"))
	data <- list(enrichedCategories=unname(rowSplit(extractSig)),
				 hasGeneSetDes=!is.null(geneSetDes), methodIsGsea=enrichMethod=="GSEA",
				 methodIsOra=enrichMethod=="ORA", fdrMethod=fdrMethod
				 )
	return(whisker.render(template, data, partials=list(statDes=statDes)))
}
