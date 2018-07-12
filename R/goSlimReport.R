goSlimReport <- function(timeStamp){
	goSlimPicPath <- paste("goslim_summary_",timeStamp,".png",sep="")
	template <- readLines(system.file("inst/templates/goSlimReport.mustache", package="WebGestaltR"))
	return(whisker.render(template, list(goSlimPicPath=goSlimPicPath)))
}
