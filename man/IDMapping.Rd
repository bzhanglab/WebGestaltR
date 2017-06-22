\name{IDMapping}
\alias{IDMapping}

\title{
	ID Mapping
}

\description{
	The IDMapping function can map one ID type supported by the WebGestaltR 
	to any other ID type supported by the WebGestaltR. This function can perform 
	the ID mapping for three types of data: a gene list, 
	a gene list with the scores and a gmt file.
}

\usage{
IDMapping(organism="hsapiens",dataType="list", inputGeneFile=NULL, 
inputGene=NULL, sourceIdType, targetIdType, collapseMethod="mean", 
mappingOutput=FALSE, outputFileName="",methodType="R",
hostName="http://www.webgestalt.org/")
}

\arguments{
  	\item{organism}{
	Currently, WebGestaltR supports 12 organisms. 
	Users can use the function \code{listOrganism} to check the available organisms.
	}
	\item{dataType}{
	The IDMapping function can perform the ID mapping for three types of data: 
	\code{list} (a gene list), \code{rnk} (a ranked gene list) and \code{gmt} (a gmt file).
	}
	\item{inputGeneFile}{
	Three types of the files are supported for uploading the data: 
	a \code{txt} file for a gene list, 
	a \code{rnk} file for a gene list with scores (separated by tab),
	a \code{gmt} file (first column is category id, second one is external link of 
	the category and others are the annotated genes. all columns are separated by tab).
	}
	\item{inputGene}{
	Two types of the R objects are supported for uploading the data: 
	an R \code{vector} object for a gene list and 
	an R \code{data.frame} object for a gene list with scores.
	}
	\item{sourceIdType}{
	The ID type of the uploaded data. 
	The supported ID type of the WebGestaltR for the selected organism 
	can be found by the function \code{listIDType}.
	}
	\item{targetIdType}{
	The target ID type for ID mapping. 
	The supported ID type of the WebGestaltR for the selected organism 
	can be found by the function \code{listIDType}.
	}
	\item{collapseMethod}{
	The method to collapse the duplicate ids for a gene list with scores 
	before ID mapping. \code{mean}, \code{median}, \code{min} and \code{max} 
	represent the mean, median, minimum and maximum of scores for the duplicate ids.
	}
	\item{mappingOutput}{
	If mappingOutput is TRUE, the mapping results will be outputted to a file.
	}
	\item{outputFileName}{
	The output file name. No extension in the file name and 
	the function will add the extension based on the input data type.
	}
	\item{methodType}{
	For the large ID mapping table (e.g. dbSNP mapping table),
	Users can use \code{R} or \code{Python} function to read it.
	Sometimes \code{Python} code is faster than the \code{R}
	code.
	If users select to use \code{Python} code to read the mapping
	table, users should first install python and the module \code{pandas}
	in the computer.
	}
	\item{hostName}{
  The server URL for accessing the data. User can use \code{listArchiveURL} function 
  to get all archive version URL.
  }
}

\value{
	The IDMapping function will output an R data.frame obejct with three types of 
	structure based on the three types of the input data. If the \code{targetIdType} 
	is one of \code{entrezgene}, \code{genesymbol} and \code{genename}, the output 
	object will contain four columns for a gene list (\code{userid}, \code{genesymbol}, 
	\code{genename} and \code{entrezgene}), five columns for a gene list with scores 
	(\code{userid}, \code{genesymbol}, \code{genename}, \code{entrezgene} and \code{score}) 
	and six columns for a gmt file (\code{geneset}, \code{link}, \code{userid}, 
	\code{genesymbol}, \code{genename} and \code{entrezgene}). If the \code{targetIdType} 
	is other ID type, the data.frame object will add one more column \code{targetid}.
}

\note{
	Because the IDMapping function will read the mapping tables from the server, 
	the running time for the WebGestaltR function will be also based on the 
	internet speed. Generally, it will take around 20 seconds to perform the ID mapping.

}

\author{
	Jing Wang
}

\examples{
#geneFile<-system.file("extdata","interestingGenes.txt",package="WebGestaltR")
#idmap <- IDMapping(organism="hsapiens",dataType="list",
#inputGeneFile=geneFile, sourceIdType="genesymbol",targetIdType="entrezgene",
#mappingOutput=FALSE,outputFileName="",methodType="R",
#hostName="http://www.webgestalt.org/")
}

\keyword{ methods }
