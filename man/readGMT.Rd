\name{readGMT}
\alias{readGMT}

\title{
	Read GMT file
}

\description{
	The readGMT function can read the GMT file and transform to an R matrix object 
	containing three columns: category ID, external link of the category 
	and the annotated genes.
}

\usage{
readGMT(gmtFile)
}

\arguments{
  \item{gmtFile}{
	The GMT file with the extension \code{gmt}.
  }
}

\author{
	Jing Wang
}

\value{
	This function will return an R matrix object containing three columns: 
	category ID, external link of the category and the annotated genes.
}

\examples{
	#gmtFile <- system.file("extdata","exampleGMT.gmt",package="WebGestaltR")
	#gmtData <- readGMT(gmtFile)
}

\keyword{ methods }
