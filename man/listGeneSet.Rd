\name{listGeneSet}
\alias{listGeneSet}

\title{
	List gene sets
}

\description{
	The listGeneSet function can list all available gene sets for 
	the selected organism in the WebGestaltR.
}

\usage{
listGeneSet(organism="hsapiens",hostName="http://www.webgestalt.org/")
}

\arguments{
  \item{organism}{
	Currently, the \code{listGeneSet} function supports 12 organisms. 
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
	The listGeneSet function can list all available gene sets for the 
	selected organism in the WebGestaltR.
}

\examples{
	#geneSet <- listGeneSet(organism="hsapiens",hostName="http://www.webgestalt.org/")
}

\keyword{ methods }
