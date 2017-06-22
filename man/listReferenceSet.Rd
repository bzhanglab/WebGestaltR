\name{listReferenceSet}
\alias{listReferenceSet}

\title{
	List reference set
}

\description{
	The listReferenceSet function can list all existing reference gene sets 
	for the selected organism in the WebGestaltR.
}

\usage{
listReferenceSet(organism="hsapiens",hostName="http://www.webgestalt.org/")
}

\arguments{
  \item{organism}{
	Currently, the \code{referenceSet} function supports 12 organisms. 
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
	The listReferenceSet function can list all existing reference gene sets 
	for the selected organism in the WebGestaltR.
}

\examples{
	#referenceSet <- listReferenceSet(organism="hsapiens",
	#hostName="http://www.webgestalt.org/")
}

\keyword{ methods }
