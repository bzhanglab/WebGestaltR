detailedGeneTable <- function(outputHtmlFile,organism,enrichMethod,fdrMethod,enrichedSig_sub,geneSetDes,geneSetNet,interestingGeneMap,enrichDatabase){

	standardId <- interestingGeneMap$standardId
	#############################################
	cat('<div id="genelist_wrap">\n',file=outputHtmlFile,append=TRUE)
	cat('<div id="genelist_des">\n',file=outputHtmlFile,append=TRUE)
				 		
	statisticDescription(enrichMethod,outputHtmlFile,fdrMethod)   ##This function include the hard code
				 		
	##########Search ID##############
	geneTableSearch(enrichedSig_sub,outputHtmlFile,geneSetDes)
	 		
  ###Show the detailed gene lists########
	plotGeneTable(outputHtmlFile,organism,enrichedSig_sub,enrichMethod,interestingGeneMap,geneSetDes,geneSetNet,enrichDatabase)

	cat("</div>\n",file=outputHtmlFile,append=TRUE)				 
	cat("</div>\n",file=outputHtmlFile,append=TRUE)
}


statisticDescription <- function(enrichMethod,outputHtmlFile,fdrMethod){
	####white the description of the statistic in the HTML page
	cat('<h4>Detailed information of the enriched categories</h4>\n',file=outputHtmlFile,append=TRUE)
	cat('<div class="tooltip_templates">The statistics\n',file=outputHtmlFile,append=TRUE)
	cat('<a class="tooltip-container"><img src="http://www.webgestalt.org/images/info.png" width="13" height="13"/>\n',file=outputHtmlFile,append=TRUE)
	cat('<span class="tooltip_content"><strong>',file=outputHtmlFile,append=TRUE)
	
	if(enrichMethod=="ORA"){
		cat("<ul><li>C: the number of reference IDs in the category</li>",file=outputHtmlFile,append=TRUE)
		cat("<li>O: the number of IDs in the user uploaded list and also in the category</li>",file=outputHtmlFile,append=TRUE)
		cat("<li>E: The expected number in the category</li>",file=outputHtmlFile,append=TRUE)
		cat("<li>R: ratio of enrichment</li>",file=outputHtmlFile,append=TRUE)
		cat("<li>PValue: p value from hyergeometric test</li>",file=outputHtmlFile,append=TRUE)
		cat("<li>FDR: FDR from ",fdrMethod,"</li></ul>",file=outputHtmlFile,append=TRUE,sep="")
		cat("</strong></span></a>\n",file=outputHtmlFile,append=TRUE,sep="")
		cat(" for the enriched categories and the IDs in the user uploaded list and also in the category are listed in the table. <br/>",file=outputHtmlFile,append=TRUE)
	}
	if(enrichMethod=="GSEA"){
		cat("<ul><li>Size: the number of IDs in the user uploaded list and also in the category</li>",file=outputHtmlFile,append=TRUE)
		cat("<li>L: the number of leading edge IDs</li>",file=outputHtmlFile,append=TRUE)
		cat("<li>ES: Enrichment Score</li>",file=outputHtmlFile,append=TRUE)
		cat("<li>NES: Normalized Enrichment Score</li>",file=outputHtmlFile,append=TRUE)
		cat("<li>PValue: p value from GSEA</li>",file=outputHtmlFile,append=TRUE)
		cat("<li>FDR: FDR from GSEA</li></ul>",file=outputHtmlFile,append=TRUE)
		cat("</strong></span></a>\n",file=outputHtmlFile,append=TRUE,sep="")
		cat(" for the enriched categories and the leading edge IDs are listed in the table.<br/>",file=outputHtmlFile,append=TRUE)
	}
				 		
	cat("</div>\n",file=outputHtmlFile,append=TRUE)
}


geneTableSearch <- function(enrichedSig_sub,outputHtmlFile,geneSetDes){
	sigID <- unique(enrichedSig_sub[,1])
	cat('<div class="search">',file=outputHtmlFile,append=TRUE)
	cat('<input id="tags" type="text" name="q" placeholder="ID Search">',file=outputHtmlFile,append=TRUE)
	cat('<ul id="searchData">',file=outputHtmlFile,append=TRUE)
	for(i in c(1:nrow(enrichedSig_sub))){
		if(!is.null(geneSetDes)){
			cat('<li><a id="aid',i-1,'" href="#',enrichedSig_sub[i,"geneset"],'">',enrichedSig_sub[i,"geneset"],' ',enrichedSig_sub[i,"description"],'</a></li>',file=outputHtmlFile,append=TRUE,sep="")
		}else{
			cat('<li><a id="aid',i-1,'" href="#',enrichedSig_sub[i,"geneset"],'">',enrichedSig_sub[i,"geneset"],'</a></li>',file=outputHtmlFile,append=TRUE,sep="")
		}
	 }
	 cat("</ul></div>\n",file=outputHtmlFile,append=TRUE,sep="")
	 cat('</div>\n',file=outputHtmlFile,append=TRUE)			 		
	 cat("<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
}


plotGeneTable <- function(outputHtmlFile,organism,enrichedSig_sub,enrichMethod,interestingGeneMap,geneSetDes,geneSetNet,enrichDatabase){
	
	standardId <- interestingGeneMap$standardId
	cat('<div id="genelist">\n',file=outputHtmlFile,append=TRUE)
	for(i in c(1:nrow(enrichedSig_sub))){
				  
		id = paste("Table",i,sep="")
		
		##Hard Code#####		  		
		if(enrichMethod=="ORA"){
			g <- enrichedSig_sub[i,"overlapID"]
		}
		
		if(enrichMethod=="GSEA"){
			g <- enrichedSig_sub[i,"leadingEdgeID"]
		}
							
		gNa <- 1
		if(!is.na(g)){  ##some enriched terms from GSEA may not have leading edge genes
			g <- unlist(strsplit(g,";"))
			x <- interestingGeneMap$mapped
			g <- x[x[,standardId] %in% g,,drop=FALSE]
			gNa <- 0
		}
		
		cat('<div><nav class="cl-effect-1">',file=outputHtmlFile,append=TRUE,sep="")
				  		
		#####Download table#######
		cat("<a style=\"cursor:pointer;color:red;\" onclick=\"$('#",id,"').tableExport({type:'csv',escape:'false'});\">Download Table</a>",file=outputHtmlFile,append=TRUE,sep="")
		
		if(!is.null(geneSetNet)){
		####Visualize Network######HARD CODE########
			if(enrichMethod=="ORA"){
				cat("<a style=\"cursor:pointer;color:blue;\" onclick=\"cytoScapeWeb('",organism,"','",enrichedSig_sub[i,1],"','",paste(g[,"genesymbol"],collapse=";"),"','');\">Network View</a>",file=outputHtmlFile,append=TRUE,sep="")
			}
			if(enrichMethod=="GSEA"){
				cat("<a style=\"cursor:pointer;color:blue;\" onclick=\"cytoScapeWeb('",organism,"','",enrichedSig_sub[i,1],"','",paste(g[,"genesymbol"],collapse=";"),"','",paste(g[,"score"],collapse=";"),"');\">Network View</a>",file=outputHtmlFile,append=TRUE,sep="")
			}
		}
				  		
		cat('</nav></div>',file=outputHtmlFile,append=TRUE,sep="")
		##################################

		###Plot table####
		cat('<table id="',id,'" class="zebra"><thead><tr>\n',file=outputHtmlFile,append=TRUE,sep="")
		cat('<th colspan="',ncol(interestingGeneMap$mapped),'">\n',file=outputHtmlFile,append=TRUE,sep="")
		
		####Hard code#####
		if(enrichMethod=="ORA"){
				elink <- linkModification(enrichDatabase,enrichedSig_sub[i,"link"],enrichedSig_sub[i,"overlapID"],interestingGeneMap)
		}
		if(enrichMethod=="GSEA"){
				elink <- linkModification(enrichDatabase,enrichedSig_sub[i,"link"],enrichedSig_sub[i,"leadingEdgeID"],interestingGeneMap)
		}
		##################
					
		if(!is.null(geneSetDes)){###Have description file
				cat('<a name="',enrichedSig_sub[i,"geneset"],'"><font color="black"><a href ="',elink,'" target="_blank"><b>FuncSet:',enrichedSig_sub[i,"geneset"],'&nbsp&nbsp&nbsp&nbsp&nbsp&nbspName:',enrichedSig_sub[i,"description"],'</b></a></font></a>\n',file=outputHtmlFile,append=TRUE,sep="")
		}else{##No description file
				cat('<a name="',enrichedSig_sub[i,"geneset"],'"><font color="black"><a href ="',elink,'" target="_blank"><b>FuncSet:',enrichedSig_sub[i,"geneset"],'</b></a></font></a>\n',file=outputHtmlFile,append=TRUE,sep="")	
		}
		cat('</th></tr>\n',file=outputHtmlFile,append=TRUE)
		cat('<tr><th colspan="',ncol(interestingGeneMap$mapped),'">\n',file=outputHtmlFile,append=TRUE,sep="")
		
		####Hard code#######	
		if(enrichMethod=="ORA"){
			cat('C=',enrichedSig_sub[i,"C"],'; O=',enrichedSig_sub[i,"O"],'; E=',round(enrichedSig_sub[i,"E"],digits=2),'; R=',round(enrichedSig_sub[i,"R"],digits=2),'; PValue=',format(enrichedSig_sub[i,"PValue"],scientific=TRUE,digits=3),'; FDR=',format(enrichedSig_sub[i,"FDR"],scientific=TRUE,digits=3),'</th></tr>\n',file=outputHtmlFile,append=TRUE,sep="")
		}
				
		if(enrichMethod=="GSEA"){
			cat('Size=',enrichedSig_sub[i,"Size"],'; L=',enrichedSig_sub[i,"leadingEdgeNum"],'; ES=',round(enrichedSig_sub[i,"ES"],digits=2),'; NES=',round(enrichedSig_sub[i,"NES"],digits=2),'; Pvalue=',format(enrichedSig_sub[i,"PValue"],scientific=TRUE,digits=3),'; FDR=',format(enrichedSig_sub[i,"FDR"],scientific=TRUE,digits=3),'</th></tr>\n',file=outputHtmlFile,append=TRUE,sep="")
		}
		###############################
		
		
		tableName <- colnames(interestingGeneMap$mapped)
		cat('<tr>\n',file=outputHtmlFile,append=TRUE)
		for(j in c(1:c(length(tableName)-1))){
			cat('<th>',tableName[j],'</th>\n',file=outputHtmlFile,append=TRUE,sep="")
		}
		cat('</tr></thead>\n',file=outputHtmlFile,append=TRUE,sep="")
		cat('<tbody>\n',file=outputHtmlFile,append=TRUE,sep="")
		if(gNa==0){ ####GSEA method may have significant categories without leading edge genes. This may also happen when users use 'TOP' method to select the enriched categories. If all categories do not annotate to any interesting genes, all selcted 'TOP' categories should only have 'NA' genes
					
			for(j in c(1:nrow(g))){
				cat('<tr>\n',file=outputHtmlFile,append=TRUE)
				a <- g[j,]
				for(k in c(1:(length(a)-1))){
					if(names(a[k])==standardId){
				 		cat('<td><a target="',standardId,'" href="',as.vector(as.matrix(a["glink"])),'">',as.vector(as.matrix(a[k])),'</a></td>\n',file=outputHtmlFile,append=TRUE,sep="")
					}else{
				 		cat('<td>',as.vector(as.matrix(a[k])),'</td>\n',file=outputHtmlFile,append=TRUE,sep="")
					}
				}
				cat('</tr>\n',file=outputHtmlFile,append=TRUE)
			}
		}
		cat('</tbody></table><br/><br/>\n',file=outputHtmlFile,append=TRUE)
	}
}