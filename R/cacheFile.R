urlToFile <- function(dataUrl) {
  result<-gsub("http://","",dataUrl)
  result<-gsub("https://","", result)
  result<-gsub("\\?[^?]+?=","_",result)
  result<-gsub("&[^&]+?=","_",result)
  result<-gsub("[:/.]","_",result)
	result<-gsub("_+", "_",result)
	return(result)
}

cacheUrl <- function(dataUrl, query=NA){
  cacheDir<-".webgestalt_cache"
  dir.create(cacheDir, showWarnings = FALSE)
  if(!is.na(query)){
    localFilePrefix<-urlToFile(paste0(dataUrl, "_", paste0(query, collapse="_")))
  }else{
    localFilePrefix<-urlToFile(dataUrl)
  }
  localFile<-paste0(cacheDir, "/", localFilePrefix, ".rds")
  if(!file.exists(localFile)){
    cat("Reading from ", dataUrl, "\n")
    if(!is.na(query)){
      response <- GET(dataUrl, query=query)
    }else{
      response <- GET(dataUrl)
    }
    
    if (response$status_code != 200) {
      return(response)
    }
    saveRDS(response, localFile)
  }else{
    cat("Reading from ", localFile, "\n")
    response <-readRDS(localFile)
  }
  
  return(response)
}
