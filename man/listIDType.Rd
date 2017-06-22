\name{listIDType}
\alias{listIDType}

\title{
	List ID type
}

\description{
	The listIDType function can list all available ID types 
	for the selected organism in the WebGestaltR.
}

\usage{
listIDType(organism="hsapiens",hostName="http://www.webgestalt.org/")
}

\arguments{
  \item{organism}{
	Currently, the \code{listIDType} function supports 12 organisms. 
	Users can use the function \code{listOrganism} to check the available organisms.
  }
  \item{hostName}{
  The server URL for accessing the data. User can use \code{listArchiveURL} function 
  to get all archive version URL.
  }
}

\author{
	Jing Wang
}

\value{
	The listIDType function can list all available ID types 
	for the selected organism in the WebGestaltR.
}

\examples{
	#idType <- listIDType(organism="hsapiens",hostName="http://www.webgestalt.org/")
}

\keyword{ methods }
