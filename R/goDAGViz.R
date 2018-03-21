goDagViz <- function(enrichedSigSub,enrichMethod,geneSetDes,geneSetDag,outputHtmlFile,outputDirectory,timeStamp,dagColor,hostName){
	####Create DAG structure using GOView########
	####Create Json file#############
	jF <- createJsonFile(enrichedSigSub,enrichMethod,geneSetDes,geneSetDag,reportDir=outputHtmlFile,outputDirectory=outputDirectory,timeStamp=timeStamp,dagColor=dagColor)

	####This is the style for saving the DAG to the file. Because google chrome does not allow to
	####read the cssRule from the local css by the javascript, we need to input this style to the javascript manually
	style <- read.table(file.path(hostName,"css","xtrace_style.txt"),header=FALSE,sep="\t",stringsAsFactors=FALSE)
	style <- as.vector(as.matrix(style))
	style <- paste(style,collapse="WJ")
	style <- paste("WJ",style,"WJ",sep="")

	re <- list(jF=jF,style=style)
	return(re)
}

createJsonFile <- function(enrichedGo,enrichMethod,desFile,dagFile,reportDir,outputDirectory,timeStamp,dagColor){
	sigGo <- unique(enrichedGo[,1])
	allGo <- unique(enrichedGo[,1])
	jsonF <- paste('[{"id":"Project_',timeStamp,'_.json","reports":[',sep="")

	if(dagColor=="continuous"){
		###HARD CODE#######
		if(enrichMethod=="ORA"){
			minF <- min(enrichedGo[,"FDR"])
			minFT <- ifelse(minF==0,-log10(2.2e-16),-log10(minF))
			colF <- colorRampPalette(c("white","red"))(128)
			myBreak <- seq(0,minFT+0.01,length.out=129)
		}

		if(enrichMethod=="GSEA"){
			x <- enrichedGo[,c("NES","FDR")]
			x[x[,2]==0,2] <- 2.2e-16
			x[,2] <- (sign(x[,1])*(-log10(x[,2])))
			minFT <- min(x[,2])
			maxFT <- max(x[,2])
			if(minFT>0){
				colF <- colorRampPalette(c("white","red"))(128)
				myBreak <- seq(0,maxFT+0.01,length.out=129)
			}else{
				if(maxFT<0){
					colF <- colorRampPalette(c("blue","white"))(128)
					myBreak <- seq(minFT-0.01,0,length.out=129)
				}else{
					if(abs(minFT)>maxFT){
						colF <- colorRampPalette(c("blue","white","red"))(256)
						myBreak <- c(seq(minFT-0.01,-0.01,length.out=128),0,seq(0.01,max(minFT)+0.01,length.out=128))
					}else{
						colF <- colorRampPalette(c("blue","white","red"))(256)
						myBreak <- c(seq(-maxFT-0.01,-0.01,length.out=128),0,seq(0.01,maxFT+0.01,length.out=128))
					}
				}
			}
		}
	}

	while(length(allGo)>0){
		g <- allGo[1]
		gdes <- desFile[desFile[,1]==g,2]
		parents <- dagFile[dagFile[,2]==g,1]

		jsonF <- paste(jsonF,'{"gID":["',g,'"],"Edge":[',sep="")
		for(i in c(1:length(parents))){
			jsonF <- paste(jsonF,'"',parents[i],'",',sep="")
		}
		jsonF <- substring(jsonF,1,nchar(jsonF)-1)
		# description may contains "'" which should be escaped
		gdes <- gsub("'", "\\\\'", gdes)
		jsonF <- paste(jsonF,'],"Name":["',gdes,'"],"Sig":["',sep="")

		if(length(which(sigGo==g))>0){
			l <- paste(reportDir,"#",g,sep="")

			####HARD CODE#######
			if(enrichMethod=="ORA"){
				gnum <- enrichedGo[enrichedGo[,1]==g,"O"]
				f <- enrichedGo[enrichedGo[,1]==g,"FDR"]
				adjp <- format(f,scientific=TRUE,digits=3)
				if(dagColor=="binary"){
					sig <- "red"
				}else{
					x <- ifelse(f==0,-log10(2.2e-16),-log10(f))
					sig <- colF[max(which(myBreak<x))]
				}
			}
			if(enrichMethod=="GSEA"){
				gnum <- enrichedGo[enrichedGo[,1]==g,"leadingEdgeNum"]
				f <- enrichedGo[enrichedGo[,1]==g,"FDR"]
				nes <- enrichedGo[enrichedGo[,1]==g,"NES"]
				adjp <- format(f,scientific=TRUE,digits=3)
				if(dagColor=="binary"){
					sig <- ifelse(nes>0,"red","blue")
				}else{
					x <- ifelse(f==0,sign(nes)*(-log10(2.2e-16)),sign(nes)*(-log10(f)))
					sig <- colF[max(which(myBreak<x))]
				}
			}
		}else{
			sig <- "white"
			l <- ""
			gnum <- ""
			adjp <- ""
		}

		jsonF <- paste(jsonF,sig,'"],"Link":["',l,'"],"gNum":["',gnum,'"],"FDR":["',adjp,'"]},',sep="")
		allGo <- allGo[-1]
		allGo <- unique(c(allGo,parents))
	}
	jsonF <- substring(jsonF,1,nchar(jsonF)-1)
	jsonF <- paste(jsonF,"]}]",sep="")
	outputFileName <- file.path(outputDirectory,paste("Project_",timeStamp,sep=""),paste("Project_",timeStamp,"_.json",sep=""))
	write.table(jsonF,file=outputFileName,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
	return(jsonF)
}
