listArchiveUrl <- function(){
	archiveUrl <- fread(input="http://www.webgestalt.org/archiveURL.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses="character",data.table=FALSE,showProgress=FALSE)
	return(archiveUrl)
}

listArchiveURL <- function(...) {
	cat("WARNING: Function listArchiveURL is deprecated and changed to listArchiveUrl!\n")
	return(listArchiveUrl(...))
}
