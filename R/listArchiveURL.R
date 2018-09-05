listArchiveUrl <- function(){
	archiveUrl <- read_tsv("http://www.webgestalt.org/archiveURL.txt", col_names=False)
	return(archiveUrl)
}

listArchiveURL <- function(...) {
	cat("WARNING: Function listArchiveURL is deprecated and changed to listArchiveUrl!\n")
	return(listArchiveUrl(...))
}
