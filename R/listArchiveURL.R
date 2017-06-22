listArchiveURL <- function(){
	archiveURL <- fread(input="http://www.webgestalt.org/archiveURL.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses="character",data.table=FALSE,showProgress=FALSE)
	return(archiveURL)
}
