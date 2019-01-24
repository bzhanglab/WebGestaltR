cacheFile <- function(hostName, paths){
  cacheDir<-".webgestalt_cache"
  dir.create(cacheDir, showWarnings = FALSE)
  localFile<-paste0(cacheDir, "/", paste0(paste(paths, collapse="_"), ".json"))
  if(!file.exists(localFile)){
    response <- GET(file.path(hostName, "api", paste(paths, collapse="/")))

    if (response$status_code != 200) {
      return(list(Succeed=FALSE, Error=webRequestError(response), jsonData=NA ))
    }
    jsonData <- content(response)
    jsonStr <-toJSON(jsonData)
    fileConn<-file(localFile)
    writeLines(jsonStr, fileConn)
    close(fileConn)
  }else{
    jsonData <-fromJSON(file=localFile)
  }
  
  return(list(Succeed=TRUE, Error=NA, jsonData=jsonData))
}

cacheFileTxt <- function(hostName, paths, query) {
  cacheDir<-".webgestalt_cache"
  dir.create(cacheDir, showWarnings = FALSE)
  localFile<-paste0(cacheDir, "/", paste0(paste(c(paths, query), collapse="_"), ".txt"))
  if(!file.exists(localFile)){
    response <- GET(file.path(hostName, "api", paste(paths, collapse="/")), query=query)
    
    if (response$status_code != 200) {
      return(list(Succeed=FALSE, Error=webRequestError(response), txtData=NA ))
    }
    txtData <- content(response)
    fileConn<-file(localFile)
    writeLines(txtData, fileConn)
    close(fileConn)
  }else{
    txtData <-readLines(localFile)
  }
  
  return(list(Succeed=TRUE, Error=NA, txtData=txtData))
}
