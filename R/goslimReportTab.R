goslimReportTab <- function(outputHtmlFile,timeStamp){
	cat('<div id="referenceGene">\n',file=outputHtmlFile,append=TRUE)
	cat('<h4>GOSlim summary for the user uploaded IDs</h4>\n',file=outputHtmlFile,append=TRUE)
	cat('Each Biological Process, Cellular Component and Molecular Function category is represented by a red, blue and green bar, repectively.<br/>\n',file=outputHtmlFile,append=TRUE)
	cat('The height of the bar represents the number of IDs in the user list and also in the category.\n',file=outputHtmlFile,append=TRUE)
	goslim <- paste("goslim_summary_",timeStamp,".png",sep="")
	cat('<img src="',goslim,'" width="100%" height="100%"/>',file=outputHtmlFile,append=TRUE,sep="")
	#cat('<iframe src="',goslim,'" width="100%" height="100%">This brower does not support PDFs. Please download the PDF to view it: <a href="',goslim,'">Download PDF</a></iframe>\n',file=outputHtmlFile,append=TRUE,sep="")
	cat("</div>\n",file=outputHtmlFile,append=TRUE)
}
