\name{formatCheck}
\alias{formatCheck}

\title{
	Check the format of the uploaded gene list or ranked gene list file or object.
}

\description{
	The formatCheck function can check the format of the gene list or 
	ranked gene list file or object uploaded to the WebGestaltR for the analysis.
}

\usage{
formatCheck(dataType="list", inputGeneFile=NULL,inputGene=NULL)
}

\arguments{
  \item{dataType}{
	Currently, the \code{formatCheck} function supports 2 data type: 
	1. "list" means the uploaded file or data is a gene list
	2. "rnk" means the uploaded file or data is a ranked gene list 
	with two columns (genes and scores).
  }
   \item{inputGeneFile}{
	The uploaded data file. 
	If the dataType is "list", the file extension should be "txt" 
	and have only one column. 
	If the dataType is "rnk", the file extension should be "rnk" 
	and have two columns (genes and scores).
  }
   \item{inputGene}{
	The uploaded R object. 
	If the dataType is "list", the R object should be a vector. 
	If the dataType is "rnk", the R object should be a data.frame 
	and have two columns (genes and scores).
  } 
}

\author{
	Jing Wang
}

\value{
	The formatCheck function will return error if the data format is incorrect. 
	Otherwise, it will return the processed data for the analysis.
}

\examples{
#geneFile<-system.file("extdata","interestingGenes.txt",package="WebGestaltR")
#interestGene<-formatCheck(dataType="list",inputGeneFile=geneFile,inputGene=NULL)
}

\keyword{ methods }
