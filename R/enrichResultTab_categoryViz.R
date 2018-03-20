enrichResultTab_categoryViz <- function(outputHtmlFile,organism,enrichMethod,fdrMethod,enrichedSig,dNum,geneSetDAG,geneSetDes,geneSetNet,outputDirectory,timeStamp,dagColor,hostName,interestingGeneMap,enrichDatabase){
	enrichedSig_sub <- extractSubSig(enrichMethod,enrichedSig,dNum)  ###extract dNum significant categories. extractSubSig contains hard code

	data <- list(hasGeneSetDag=!is.null(geneSetDAG), methodIsGsea=enrichMethod=="GSEA",
				dagColocrIsBinary=dagColor=="binary", standardId=interestingGeneMap$standardId,
				hasGeneSetDes=!is.null(geneSetDes), methodIsOra=enrichMethod=="ORA"
				)

	dagInfo <- NULL
	enrichedCategories <- NULL
	if(!is.null(geneSetDAG)){
		dagInfo <- goDAGViz(enrichedSig_sub,enrichMethod,geneSetDes,geneSetDAG,outputHtmlFile,outputDirectory,timeStamp,dagColor,hostName)
		if(!is.null(enrichedSig) && !is.null(geneSetDAG)){
			data$jId <- paste("Project_",timeStamp,"_.json",sep="")
			data$jF <- dagInfo$jF
			data$style <- dagInfo$style
		}
		data$enrichedCategories <- enrichedCategories
	}else{
		#########Create table to summary enriched results###########
		if(enrichMethod=="ORA"){
			enrichedCategories <- cbind(enrichedSig_sub, "color"=rep("black", nrow(enrichedSig_sub)))
		} else if(enrichMethod=="GSEA"){
			color <- vector("character", nrow(enrichedSig_sub))
			color[enrichedSig_sub[,"NES"]>0] <- "red"
			color[enrichedSig_sub[,"NES"]<=0] <- "blue"
			enrichedCategories <- cbind(enrichedSig_sub, "color"=color)
		}
		enrichedCategories[,"FDR"] <- format(enrichedCategories[, "FDR"], scientific=TRUE, digits=3)
		data$enrichedCategories <- unname(rowSplit(enrichedCategories))
	}

	data$detailedGeneTableContent <- detailedGeneTable(outputHtmlFile,organism,enrichMethod,fdrMethod,enrichedSig_sub,geneSetDes,geneSetNet,interestingGeneMap,enrichDatabase)

	template <- readLines(system.file("inst/templates/enrichResultTab.mustache", package="WebGestaltR"))
	return(whisker.render(template, data))
}

extractSubSig <- function(enrichMethod,enrichedSig,dNum){
	###extract the dNum siganificant categories from the enrichedSig
	###Hard Code######
	if(enrichMethod=="ORA"){
		enrichedSig <- enrichedSig[order(enrichedSig[,"FDR"],enrichedSig[,"PValue"]),]
		if(nrow(enrichedSig)>dNum){
			enrichedSig_sub <- enrichedSig[1:dNum,]
		}else{
			enrichedSig_sub <- enrichedSig
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

		enrichedSig_sub <- rbind(x,y)
	}
	return(enrichedSig_sub)
}


enrichResult_Others <- function(enrichMethod,enrichedSig,geneSetDes,fdrMethod,dNum){
	enrichedSig_sub <- extractSubSig(enrichMethod,enrichedSig,dNum)
	extractSig <- data.frame(id=enrichedSig_sub[,"geneset"],link=enrichedSig_sub[,"link"],stringsAsFactors=FALSE)

	if(enrichMethod=="ORA"){
		extractSig$statistic <- paste("C=",enrichedSig_sub[,"C"],"; O=",enrichedSig_sub[,"O"],"; E=",round(enrichedSig_sub[,"E"],digits=2),"; R=",round(enrichedSig_sub[,"R"],digits=2),"; PValue=",format(enrichedSig_sub[,"P Value"],scientific=TRUE,digits=3),"; FDR=",format(enrichedSig_sub[,"FDR"],scientific=TRUE,digits=3),sep="")
		extractSig$genes <- gsub(",|;"," ",enrichedSig_sub[,"overlapID"])
		extractSig$color <- rep("black", nrow(extractSig))
	} else if (enrichMethod=="GSEA"){
		extractSig$statistic <- paste("Size=",enrichedSig_sub[,"Size"],"; L=",enrichedSig_sub[,"leadingEdgeNum"],"; ES=",round(enrichedSig_sub[,"ES"],digits=2),"; NES=",round(enrichedSig_sub[,"NES"],digits=2),"; P value=",format(enrichedSig_sub[,"PValue"],scientific=TRUE,digits=3),"; FDR=",format(enrichedSig_sub[,"FDR"],scientific=TRUE,digits=3),sep="")
		extractSig$genes <- gsub(",|;"," ",enrichedSig_sub[,"leadingEdgeID"])
		extractSig$color <- sapply(enrichedSig_sub[,"NES"], function(x) {if(x>0) {return("red")} else {return("blue")}})
	}

	if(!is.null(geneSetDes)){
		extractSig$name <- enrichedSig_sub[,"description"]
	}

	template <- readLines(system.file("inst/templates/enrichResultOthers.mustache", package="WebGestaltR"))
	statDes <- readLines(system.file("inst/templates/enrichResultStat.mustache", package="WebGestaltR"))
	data <- list(enrichedCategories=unname(rowSplit(extractSig)),
				 hasGeneSetDes=!is.null(geneSetDes), methodIsGsea=enrichMethod=="GSEA",
				 methodIsOra=enrichMethod=="ORA", fdrMethod=fdrMethod
				 )
	return(whisker.render(template, data, partials=list(statDes=statDes)))
}
