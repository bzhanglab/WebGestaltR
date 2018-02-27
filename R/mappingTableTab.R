mappingTableTab <- function(outputHtmlFile,interestingGeneMap){
	standardId <- interestingGeneMap$standardId
	cat('<div id="interestGene">\n',file=outputHtmlFile,append=TRUE)
	cat('<div id="mappedGene">\n',file=outputHtmlFile,append=TRUE)
	cat('<h4>Mapped User IDs</h4>\n',file=outputHtmlFile,append=TRUE)
	tableName <- colnames(interestingGeneMap$mapped)
	cat('<table class="zebra"><thead><tr>\n',file=outputHtmlFile,append=TRUE)
	####plot the table name##########
	for(i in c(1:(length(tableName)-1))){
		cat('<th>',tableName[i],'</th>\n',file=outputHtmlFile,append=TRUE,sep="")
	}
	cat('</tr></thead>\n',file=outputHtmlFile,append=TRUE,sep="")

	cat('<tbody>\n',file=outputHtmlFile,append=TRUE,sep="")
	for(i in c(1:nrow(interestingGeneMap$mapped))){
		cat('<tr>\n',file=outputHtmlFile,append=TRUE)
		a <- interestingGeneMap$mapped[i,]
		for(j in c(1:(length(a)-1))){
			if(names(a[j])==standardId){
				cat('<td><a target="',standardId,'" href="',as.vector(as.matrix(a["glink"])),'">',as.vector(as.matrix(a[j])),'</a></td>\n',file=outputHtmlFile,append=TRUE,sep="")
			}else{
				cat('<td>',as.vector(as.matrix(a[j])),'</td>\n',file=outputHtmlFile,append=TRUE,sep="")
			}
		}
		cat('</tr>\n',file=outputHtmlFile,append=TRUE)
	}

	cat('</tbody></table>\n',file=outputHtmlFile,append=TRUE)
	cat("</div>\n",file=outputHtmlFile,append=TRUE)###For mappedGene

	cat('<div id="empty">&nbsp&nbsp&nbsp&nbsp</div>\n',file=outputHtmlFile,append=TRUE)

	cat('<div id="unmappedGene">\n',file=outputHtmlFile,append=TRUE)
	cat("<h4>User IDs not mapped</h4>\n",file=outputHtmlFile,append=TRUE)

	cat('<table class="zebra"><thead><tr>\n',file=outputHtmlFile,append=TRUE)

	cat('<th>userid</th>\n',file=outputHtmlFile,append=TRUE)

	cat('</tr></thead>\n',file=outputHtmlFile,append=TRUE,sep="")
	cat('<tbody>\n',file=outputHtmlFile,append=TRUE,sep="")

	for(i in c(1:length(interestingGeneMap$unmapped))){
		cat('<tr>\n',file=outputHtmlFile,append=TRUE)
		cat('<td>',interestingGeneMap$unmapped[i],'</td>\n',file=outputHtmlFile,append=TRUE,sep="")
		cat('</tr>\n',file=outputHtmlFile,append=TRUE)
	}

	cat('</tbody></table>\n',file=outputHtmlFile,append=TRUE)
	cat("</div>\n",file=outputHtmlFile,append=TRUE)  ##For empty part
	cat("</div>\n",file=outputHtmlFile,append=TRUE)  ##For interestGene part
}
