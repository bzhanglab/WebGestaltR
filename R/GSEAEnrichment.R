GSEAEnrichment <- function(hostName,outputDirectory,projectName,geneRankList,geneSet,collapseMethod="mean",minNum=10,maxNum=500,sigMethod="fdr",fdrThr=0.05,topThr=10,perNum=1000,lNum=20,is.output=TRUE,keepGSEAFolder=FALSE){
	GSEAJarFile <- list.files(path=outputDirectory,pattern="gsea.*\\.jar$")
	if(length(GSEAJarFile)==0){
		cat("No GSEA java jar file can be found in the current working directory. The function will download the GSEA java jar file to the ",outputDirectory,". The copyright of the GSEA java jar file belongs to the broad institute (http://software.broadinstitute.org/gsea/index.jsp).\n",sep="")
		GSEAJarFile <- file.path(outputDirectory,"gsea.jar")
		if(length(grep("^file://", hostName, perl=TRUE))==1){
			# copy from local file instead of downloading from remote URL
			local.file <- gsub("^file://", "", hostName, perl=TRUE)
			file.copy(file.path(local.file, "gsea.jar"), GSEAJarFile)
		} else {
			download.file(file.path(hostName,"gsea.jar"),GSEAJarFile,mode="wb")
		}
	}else{
		GSEAJarFile <- file.path(outputDirectory,GSEAJarFile)
	}

	if(length(grep(" ",GSEAJarFile))>0){
		GSEAJarFile <- gsub(" ","\\ ",GSEAJarFile,fixed=TRUE)
	}

	projectFolder <- file.path(outputDirectory,paste("Project_",projectName,sep=""))
	if(!dir.exists(projectFolder)){
		dir.create(projectFolder)
	}

	geneRankList[,1] <- as.character(geneRankList[,1])
	geneSet[,3] <- as.character(geneSet[,3])

	x <- geneSet[geneSet[,3] %in% geneRankList[,1],,drop=FALSE]

	geneSetNum <- tapply(x[,3],x[,1],length)
	geneSetNum <- geneSetNum[geneSetNum>=minNum & geneSetNum<=maxNum]
	if(length(geneSetNum)==0){
		error <- paste("ERROR: The number of annotated IDs for all functional categories are not from ",minNum," to ",maxNum," for the GSEA enrichment method.",sep="")
		cat(error)
		return(error)
	}

	a <- tapply(geneRankList[,2],geneRankList[,1],collapseMethod,na.rm=TRUE)
	geneRankList <- data.frame(geneid=names(a),score=a,stringsAsFactors=FALSE)

	gsea_rnk <- file.path(projectFolder,paste("Project_",projectName,"_GSEA.rnk",sep=""))
	write.table(geneRankList,file=gsea_rnk,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)

	gsea_gmt <- file.path(projectFolder,paste("Project_",projectName,"_GSEA.gmt",sep=""))
	x <- tapply(geneSet[,3],geneSet[,1],paste,collapse="\t")
	y <- unique(geneSet[,c(1,2)])
	x <- x[y[,1]]
	g <- data.frame(geneset=y[,1],link=y[,2],gene=x,stringsAsFactors=FALSE)
	write.table(g,file=gsea_gmt,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)

	outputF <- file.path(projectFolder,paste("Project_",projectName,"_GSEA/",sep=""))
	if(!dir.exists(outputF)){
		dir.create(outputF)
	}

	comm <- paste("java -Xmx4096m -cp ",GSEAJarFile," xtools.gsea.GseaPreranked -gmx ",gsea_gmt," -collapse false -mode Max_probe -norm meandiv -nperm ",perNum," -rnk ",gsea_rnk," -scoring_scheme weighted -rpt_label Project_",projectName," -include_only_symbols true -make_sets true -plot_top_x ",lNum," -rnd_seed timestamp -set_max ",maxNum," -set_min ",minNum," -zip_report false -out ",outputF," -gui false",sep="")

	eI <- system(comm,ignore.stderr=TRUE)
	if(eI>0){
		error <- "ERROR: please install the Java or set the correct path of Java source file in the system."
		cat(error)
		return(error)
	}

	gseaR <- readGSEA(outputF)
	gseaR_pos <- gseaR$positiveP
	gseaR_neg <- gseaR$negativeP
	if(is.null(gseaR_pos) && is.null(gseaR_neg)){
		error <- "No enriched set is outputted from GSEA"
		cat(error)
		return(error)
	}else{
		if(sigMethod=="fdr"){
			if(!is.null(gseaR_pos)){
				gseaR_posSig <- gseaR_pos[gseaR_pos[,"FDR"]<fdrThr,]
			}else{
				gseaR_posSig <- NULL
			}

			if(!is.null(gseaR_neg)){
				gseaR_negSig <- gseaR_neg[gseaR_neg[,"FDR"]<fdrThr,]
			}else{
				gseaR_negSig <- NULL
			}

			sig <- rbind(gseaR_posSig,gseaR_negSig)
			if(nrow(sig)==0){
				cat("No significant set is identified based on FDR ",fdrThr,"!",sep="")
				removeFolder(projectFolder,is.output=is.output,keepGSEAFolder=keepGSEAFolder)
				return(NULL)
			}else{
				sig <- mappingName(sig,geneSet)
				sig <- sig[order(sig[,"FDR"],sig[,"NES"]),]
				removeFolder(projectFolder,is.output=is.output,keepGSEAFolder=keepGSEAFolder)
				return(sig)
			}
		}else{
			if(!is.null(gseaR_pos)){
				gseaR_pos <- gseaR_pos[order(gseaR_pos[,"FDR"],gseaR_pos[,"PValue"]),]
				if(nrow(gseaR_pos)>topThr){
					gseaR_posSig <- gseaR_pos[1:topThr,]
				}else{
					gseaR_posSig <- gseaR_pos
				}
			}else{
				gseaR_posSig <- NULL
			}

			if(!is.null(gseaR_neg)){
				gseaR_neg <- gseaR_neg[order(gseaR_neg[,"FDR"],gseaR_neg[,"PValue"]),]
				if(nrow(gseaR_neg)>topThr){
					gseaR_negSig <- gseaR_neg[1:topThr,]
				}else{
					gseaR_negSig <- gseaR_neg
				}
			}else{
				gseaR_negSig <- NULL
			}

			sig <- rbind(gseaR_posSig,gseaR_negSig)
			sig <- mappingName(sig,geneSet)
			sig <- sig[order(sig[,"FDR"],sig[,"PValue"]),]

			removeFolder(projectFolder,is.output=is.output,keepGSEAFolder=keepGSEAFolder)
			return(sig)
		}
	}
}

readGSEA <- function(gseaFolder){
	positiveP <- data.frame(geneset="",link="",Size=0,ES=0,NES=0,PValue=0,FDR=0,leadingEdgeNum=0,leadingEdgeID="",stringsAsFactors=FALSE)
	rei <- 1

	negativeP <- data.frame(geneset="",link="",Size=0,ES=0,NES=0,PValue=0,FDR=0,leadingEdgeNum=0,leadingEdgeID="",stringsAsFactors=FALSE)
	nei <- 1

	subFile <- list.files(gseaFolder,pattern="GseaPreranked")
	subF <- file.path(gseaFolder,subFile)

	positiveF <- list.files(subF,"gsea_report_for_na_pos_",full.names=TRUE)
	positiveF <- positiveF[grep("xls",positiveF,fixed=TRUE)]

	negativeF <- list.files(subF,"gsea_report_for_na_neg_",full.names=TRUE)
	negativeF <- negativeF[grep("xls",negativeF,fixed=TRUE)]

	positiveD <- fread(input=positiveF,header=TRUE,sep="\t",stringsAsFactors=FALSE,data.table=FALSE,showProgress=FALSE)

	leadingF <- list.files(subF,pattern="xls")

	##if no gene is annotated to any category, GSEA will not output any category
	if(nrow(positiveD)>0){
		for(j in c(1:nrow(positiveD))){
			positiveP[j,1] <- positiveD[j,1]
			positiveP[j,2] <- positiveD[j,2]
			positiveP[j,3] <- positiveD[j,4]
			positiveP[j,4] <- positiveD[j,5]
			positiveP[j,5] <- positiveD[j,6]
			positiveP[j,6] <- positiveD[j,7]
			positiveP[j,7] <- positiveD[j,8]
			lF <- paste(positiveD[j,1],".xls",sep="")
			if(length(which(leadingF==lF))>0){
				x <- readLeadingEdge(subF,positiveD[j,1])
				positiveP[j,8] <- x$geneListNum
				positiveP[j,9] <- x$genelist
			}else{
				positiveP[j,8] <- NA
				positiveP[j,9] <- NA
			}
		}
		positiveP <- positiveP[,c(1,2,3,8,4:7,9)]
		positiveP <- positiveP[!is.na(positiveP[,"NES"]),]     ###GSEA may generate some terms with NA NES and Pvalue
		if(nrow(positiveP)==0){
			positiveP <- NULL
		}
	}else{
		positiveP <- NULL
	}

	negativeD <- fread(input=negativeF,header=TRUE,sep="\t",stringsAsFactors=FALSE,data.table=FALSE,showProgress=FALSE)
	if(nrow(negativeD)>0){
		for(j in c(1:nrow(negativeD))){
			negativeP[j,1] <- negativeD[j,1]
			negativeP[j,2] <- negativeD[j,2]
			negativeP[j,3] <- negativeD[j,4]
			negativeP[j,4] <- negativeD[j,5]
			negativeP[j,5] <- negativeD[j,6]
			negativeP[j,6] <- negativeD[j,7]
			negativeP[j,7] <- negativeD[j,8]
			lF <- paste(negativeD[j,1],".xls",sep="")
			if(length(which(leadingF==lF))>0){
				x <- readLeadingEdge(subF,negativeD[j,1])
				negativeP[j,8] <- x$geneListNum
				negativeP[j,9] <- x$genelist
			}else{
				negativeP[j,8] <- NA
				negativeP[j,9] <- NA
			}
		}
		negativeP <- negativeP[,c(1,2,3,8,4:7,9)]
		negativeP <- negativeP[!is.na(negativeP[,"NES"]),]
		if(nrow(negativeP)==0){
			negativeP <- NULL
		}
	}else{
		negativeP <- NULL
	}

	re <- list(positiveP=positiveP,negativeP=negativeP)
	return(re)
}

readLeadingEdge <- function(dir,sigModule){
	moduleF <- fread(input=file.path(dir,paste(sigModule,".xls",sep="")),header=TRUE,sep="\t",stringsAsFactors=FALSE,data.table=FALSE,showProgress=FALSE)
	sigGene <- moduleF[moduleF[,8]=="Yes",2]
	sigGeneNum <- length(sigGene)
	sigGeneL <- paste(sigGene,collapse=";")
	re <- list(geneListNum=sigGeneNum,genelist=sigGeneL)
	return(re)
}

removeFolder <- function(dir,is.output=TRUE,keepGSEAFolder=FALSE){  ###For the GSEA enrichment analysis
	if(is.output==FALSE){
		comm <- paste("rm -rf ",dir,sep="")
		system(comm)
	}else{
		if(keepGSEAFolder==FALSE){
			all <- list.files(path=dir,full.names=TRUE)
			folder <- all[file.info(all)$isdir]
			for(i in c(1:length(folder))){
				comm <- paste("rm -rf ",folder[i],sep="")
				system(comm)
			}
		}
	}
}

mappingName <- function(gseaResult,geneSet){
####GSEA will automatically change all ID to upper-case and change link to others. We need to change ID back and add the link
	geneSetName <- unique(geneSet[,c(1,2)])
	geneSetName <- cbind(geneSetName,toupper(geneSetName[,1]))
	colnames(geneSetName) <- c("geneset1","link1","geneset")
	mE <- merge(gseaResult,geneSetName,by="geneset")
	mE <- mE[,c("geneset1","link1",colnames(gseaResult)[3:ncol(gseaResult)])]
	colnames(mE)[1:2] <- c("geneset","link")
	mE[,1] <- as.character(mE[,1])
	mE[,2] <- as.character(mE[,2])
	return(mE)
}
