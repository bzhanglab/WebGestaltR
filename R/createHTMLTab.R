createHTMLTab <- function(outputHtmlFile,enrichedSig,interestingGeneMap){

	standardId <- interestingGeneMap$standardId
	####################create the Tabs################################
	 cat('<div id ="tabs"> <ul>\n',file=outputHtmlFile,append=TRUE)
	 cat('<li><a href="#summary">Summary</a></li>\n',file=outputHtmlFile,append=TRUE)
	 cat('<li><a href="#interestGene">User ID Mapping Table</a></li>\n',file=outputHtmlFile,append=TRUE)
	 if(standardId=="entrezgene"){
	 	cat('<li><a href="#referenceGene">GOSlim Summary</a></li>\n',file=outputHtmlFile,append=TRUE)
	 }
	 if(!is.null(enrichedSig)){
		cat('<li><a href="#result">Enrichment Results</a></li>\n',file=outputHtmlFile,append=TRUE)
	 }
			
	 cat("</ul>\n",file=outputHtmlFile,append=TRUE)
}