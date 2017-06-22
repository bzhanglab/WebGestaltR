\name{listOrganism}
\alias{listOrganism}

\title{
	List organisms
}

\description{
	The listOrganism function can list all supported organisms in the WebGestaltR.
}

\usage{
listOrganism(hostName="http://www.webgestalt.org/")
}


\arguments{
  \item{hostName}{
  The server URL for accessing the data. User can use \code{listArchiveURL} function 
  to get all archive version URL.
  }
}

\author{
	Jing Wang
}

\value{
	The listOrganism function can list all supported organisms in the WebGestaltR.
}

\examples{
	#organism <- listOrganism(hostName="http://www.webgestalt.org/")
}

\keyword{ methods }
