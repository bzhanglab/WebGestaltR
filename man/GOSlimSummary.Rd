\name{GOSlimSummary}
\alias{GOSlimSummary}

\title{
	Summary the gene list based on the GO (Gene Ontology) Slim data sets
}

\description{
	This function can summary the gene list based on the biological process, 
	cellular component and molecular function ontologies of the GO Slim data sets. 
	The summary result will be plotted as three bar plots and outputted to the PDF file.
}

\usage{
GOSlimSummary(organism, genelist, outputFile, outputType, hostName)
}

\arguments{
  \item{organism}{
	Currently, GOSlimSummary supports 12 organisms. 
	Users can use the function \code{listOrganism} to check the available organisms.
  }
  \item{genelist}{
	An R \code{vector} object containing a gene list. 
	GOSlimSummary only supports NCBI EntrezGene ID for the summary. 
	For other ID types, please first use \code{IDMapping} function 
	to map to the EntrezGene ID.
  }
  \item{outputFile}{
	The output file name.
  }
  \item{outputType}{
	The output file extension that can be \code{pdf}, \code{png}, or \code{bmp}.
  }
  \item{hostName}{
  The server URL for accessing the data. User can use \code{listArchiveURL} function 
  to get all archive version URL.
  }
}

\note{
	Because the GOSlimSummary function will read the GO Slim data from the server, 
	the running time for this function will be also based on the internet speed. 
	Generally, it will take around 20 seconds to perform the summary analysis.
}

\author{
	Jing Wang
}

\value{
	GOSlimSummary function will return a high-resolution image 
	containing the summary information of the uploaded gene list 
	based on the biological process (red bar plot), cellular component (blue bar plot) 
	and molecular function (green bar plot) ontologies.
}

\examples{
#geneFile<-system.file("extdata","GOSlimExample.txt",package="WebGestaltR")
#geneList<-read.table(geneFile,header=FALSE,sep="\t",stringsAsFactors=FALSE)
#geneList<-as.vector(as.matrix(geneList))
#outputFile<-paste(getwd(),"/GOSlimSummary",sep="")
#GOSlimSummary(organism="hsapiens", genelist=geneList, 
#outputFile=outputFile, outputType="pdf")
}

\keyword{ methods }
