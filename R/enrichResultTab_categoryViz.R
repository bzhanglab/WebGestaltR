enrichResultTab_categoryViz <- function(outputHtmlFile,organism,enrichMethod,fdrMethod,enrichedSig,dNum,geneSetDAG,geneSetDes,geneSetNet,outputDirectory,timeStamp,dagColor,hostName,interestingGeneMap,enrichDatabase){	
	 standardId <- interestingGeneMap$standardId
	 enrichedSig_sub <- extractSubSig(enrichMethod,enrichedSig,dNum)  ###extract dNum significant categories. extractSubSig contains hard code
	 cat('<div id="dag">\n',file=outputHtmlFile,append=TRUE)
	 godagInfo <- NULL
	 if(!is.null(geneSetDAG)){
			godagInfo <- goDAGViz(enrichedSig_sub,enrichMethod,geneSetDes,geneSetDAG,outputHtmlFile,outputDirectory,timeStamp,dagColor,hostName)  ###GODAGViz contains HARD code for the description for the GSEA method
 	 }else{
	 #########Create table to summary enriched results###########
	 		categoryTableViz(enrichMethod,standardId,geneSetDes,enrichedSig_sub,outputHtmlFile)
	 }
	 cat("</div>\n",file=outputHtmlFile,append=TRUE)  ##For dag div
	 
	 cat('<div id="empty1">&nbsp&nbsp&nbsp&nbsp</div>\n',file=outputHtmlFile,append=TRUE)
	 detailedGeneTable(outputHtmlFile,organism,enrichMethod,fdrMethod,enrichedSig_sub,geneSetDes,geneSetNet,interestingGeneMap,enrichDatabase) 
	 cat("</div>\n",file=outputHtmlFile,append=TRUE)
	 return(godagInfo)
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



categoryTableViz <- function(enrichMethod,standardId,geneSetDes,enrichedSig_sub,outputHtmlFile){
	cat('<h4>Summary of the enriched categories</h4>\n',file=outputHtmlFile,append=TRUE)
	
	####Hard Code#######	
	if(enrichMethod=="ORA"){
			cat('This table lists the enriched categories, number of ',standardId,' in the user uploaded list and also in the categories and FDR.<br/><br/>\n',file=outputHtmlFile,append=TRUE,sep="")
	}
	if(enrichMethod=="GSEA"){
			cat('This table lists the enriched categories, number of ',standardId,' in the user uploaded list and also in the categories and FDR. Categories with red color represent the positively related categories while categories with blue color represent the negatively related categories. <br/><br/>\n',file=outputHtmlFile,append=TRUE,sep="")
	}
					
	cat('<table class="zebra"><thead><tr>\n',file=outputHtmlFile,append=TRUE)
	if(!is.null(geneSetDes)){
			cat('<th>FuncSet</th><th>Name</th><th>#',standardId,'</th><th>FDR</th>\n',file=outputHtmlFile,append=TRUE)
			if(enrichMethod=="ORA"){
				 	extractSig <- enrichedSig_sub[,c("geneset","description","O","FDR")]
			}
			if(enrichMethod=="GSEA"){
				 	extractSig <- enrichedSig_sub[,c("geneset","description","Size","FDR")]
			}
				 		
	 }else{
			cat('<th>FuncSet</th><th>#',standardId,'</th><th>FDR</th>\n',file=outputHtmlFile,append=TRUE)
			if(enrichMethod=="ORA"){
				 	extractSig <- enrichedSig_sub[,c("geneset","O","FDR")]
			}
			if(enrichMethod=="GSEA"){
				 	extractSig <- enrichedSig_sub[,c("geneset","Size","FDR")]
			}
		}
				 	
		cat('</tr></thead>\n',file=outputHtmlFile,append=TRUE,sep="")
		cat('<tbody>\n',file=outputHtmlFile,append=TRUE,sep="")
			
		for(i in c(1:nrow(extractSig))){
			cat('<tr>\n',file=outputHtmlFile,append=TRUE)
			a <- extractSig[i,]
				 				
			###Set the row color for GSEA###
			####HARD CODE####
			if(enrichMethod=="ORA"){
				 	color <- "black"
			}
			if(enrichMethod=="GSEA"){
				 	if(enrichedSig_sub[i,"NES"]>0){
				 		color <- "red"
				 	}else{
				 		color <- "blue"
				 	}
			}
				 				
			for(j in c(1:length(a))){
				 if(j==1){
				 		cat('<td><a href="#',as.vector(as.matrix(a[j])),'" style="color:',color,'">',as.vector(as.matrix(a[j])),'</a></td>\n',file=outputHtmlFile,append=TRUE,sep="")
				 }else{
				 		if(j==length(a)){
				 			cat('<td><font color="',color,'">',format(as.vector(as.matrix(a[j])),scientific=TRUE,digits=3),'</font></td>\n',file=outputHtmlFile,append=TRUE,sep="")
				 		}else{
					 		cat('<td><font color="',color,'">',as.vector(as.matrix(a[j])),'</font></td>\n',file=outputHtmlFile,append=TRUE,sep="")
					 	}
					}
			}
			cat('</tr>\n',file=outputHtmlFile,append=TRUE)
		}
		cat('</tbody></table>\n',file=outputHtmlFile,append=TRUE)
}


enrichResult_Others <- function(outputHtmlFile,enrichMethod,enrichedSig,geneSetDes,fdrMethod,dNum){
	enrichedSig_sub <- extractSubSig(enrichMethod,enrichedSig,dNum)
	statisticDescription(enrichMethod,outputHtmlFile,fdrMethod)
	cat('<div width="100%">\n',file=outputHtmlFile,append=TRUE)
	cat('<table class="zebra"><thead><tr>\n',file=outputHtmlFile,append=TRUE)
	###Hard Code######
	if(enrichMethod=="ORA"){
			if(!is.null(geneSetDes)){
					cat('<th width="10%">FuncSet</th><th width="20%">Name</th><th width="30%">Statistic</th><th width="40%">UploadedID</th></tr></thead>\n',file=outputHtmlFile,append=TRUE)
				 	extractSig <- data.frame(id=enrichedSig_sub[,"geneset"],name=enrichedSig_sub[,"description"],link=enrichedSig_sub[,"link"],statistic=paste("C=",enrichedSig_sub[,"C"],"; O=",enrichedSig_sub[,"O"],"; E=",round(enrichedSig_sub[,"E"],digits=2),"; R=",round(enrichedSig_sub[,"R"],digits=2),"; PValue=",format(enrichedSig_sub[,"PValue"],scientific=TRUE,digits=3),"; FDR=",format(enrichedSig_sub[,"FDR"],scientific=TRUE,digits=3),sep=""),genes=gsub(","," ",enrichedSig_sub[,"overlapID"]),stringsAsFactors=FALSE)
			}else{
				 	cat('<th width="10%">FuncSet</th width="30%"><th>Statistic</th><th width="60%">UploadedID</th></tr></thead>\n',file=outputHtmlFile,append=TRUE)
				 	extractSig <- data.frame(id=enrichedSig_sub[,"geneset"],link=enrichedSig_sub[,"link"],statistic=paste("C=",enrichedSig_sub[,"C"],"; O=",enrichedSig_sub[,"O"],"; E=",round(enrichedSig_sub[,"E"],digits=2),"; R=",round(enrichedSig_sub[,"R"],digits=2),"; PValue=",format(enrichedSig_sub[,"PValue"],scientific=TRUE,digits=3),"; FDR=",format(enrichedSig_sub[,"FDR"],scientific=TRUE,digits=3),sep=""),genes=gsub(","," ",enrichedSig_sub[,"overlapID"]),stringsAsFactors=FALSE)
			}
				 	
		}
		
		if(enrichMethod=="GSEA"){
				cat("Categories with red color represent positively related categories while categories with blue color represent negatively related categories.<br/>\n",file=outputHtmlFile,append=TRUE)
				if(!is.null(geneSetDes)){
				 	cat('<th width="10%">FuncSet</th><th width="20%">Name</th><th width="30%">Statistic</th><th width="40%">UploadedID</th></tr></thead>\n',file=outputHtmlFile,append=TRUE)
				 	extractSig <- data.frame(id=enrichedSig_sub[,"geneset"],name=enrichedSig_sub[,"description"],link=enrichedSig_sub[,"link"],statistic=paste("Size=",enrichedSig_sub[,"Size"],"; L=",enrichedSig_sub[,"leadingEdgeNum"],"; ES=",round(enrichedSig_sub[,"ES"],digits=2),"; NES=",round(enrichedSig_sub[,"NES"],digits=2),"; PValue=",format(enrichedSig_sub[,"PValue"],scientific=TRUE,digits=3),"; FDR=",format(enrichedSig_sub[,"FDR"],scientific=TRUE,digits=3),sep=""),genes=gsub(","," ",enrichedSig_sub[,"leadingEdgeID"]),stringsAsFactors=FALSE)
				 }else{
				 	 cat('<th width="10%">FuncSet</th width="30%"><th>Statistic</th><th width="60%">UploadedID</th></tr></thead>\n',file=outputHtmlFile,append=TRUE)
				 	 extractSig <- data.frame(id=enrichedSig_sub[,"geneset"],link=enrichedSig_sub[,"link"],statistic=paste("Size=",enrichedSig_sub[,"Size"],"; L=",enrichedSig_sub[,"leadingEdgeNum"],"; ES=",round(enrichedSig_sub[,"ES"],digits=2),"; NES=",round(enrichedSig_sub[,"NES"],digits=2),"; Pvalue=",format(enrichedSig_sub[,"PValue"],scientific=TRUE,digits=3),"; FDR=",format(enrichedSig_sub[,"FDR"],scientific=TRUE,digits=3),sep=""),genes=gsub(","," ",enrichedSig_sub[,"leadingEdgeID"]),stringsAsFactors=FALSE)
				 }
		}
		cat('<tbody>\n',file=outputHtmlFile,append=TRUE,sep="")

		for(i in c(1:nrow(enrichedSig_sub))){
				cat("<tr>\n",file=outputHtmlFile,append=TRUE)
				###HARD CODE####
				if(enrichMethod=="ORA"){
					color <- "black"
				}
				 
				if(enrichMethod=="GSEA"){
					if(enrichedSig_sub[i,"NES"]>0){
						color <- "red"
					}else{
						color <- "blue"
					}
				}
				
				id <- extractSig[i,1]
				link <- extractSig[i,"link"]
				stat <- extractSig[i,"statistic"]
				genes <- gsub(";"," ",extractSig[i,"genes"])
							
				if(!is.null(geneSetDes)){
					name <- extractSig[i,"name"]
					cat('<td><a href="',link,'" target="_blank" style="color:',color,'">',id,'</a></font></td><td><font color="',color,'">',name,'</font></td><td><font color="',color,'">',stat,'</font></td><td><font color="',color,'">',genes,'</font></td>\n',file=outputHtmlFile,append=TRUE,sep="")
				}else{
					cat('<td><a href="',link,'" target="_blank" style="color:',color,'">',id,'</a></font></td><td><font color="',color,'">',stat,'</font></td><td><font color="',color,'">',genes,'</font></td>\n',file=outputHtmlFile,append=TRUE,sep="")
				}
				cat("</tr>\n",file=outputHtmlFile,append=TRUE)
		}
		cat('</tbody></table></div><br/><br/>\n',file=outputHtmlFile,append=TRUE)
}
