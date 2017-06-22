\name{WebGestaltR_batch}
\alias{WebGestaltR_batch}

\title{
	Batch analysis for the multiple gene or ranked gene lists
}

\description{
	The WebGestaltR_batch function can perform the batch analysis for multiple 
	gene or ranked gene lists.
}

\usage{
WebGestaltR_batch(interestGeneFolder=NULL,interestGeneType=NULL,enrichMethod="ORA",
organism="hsapiens",enrichDatabase="geneontology_Biological_Process",
enrichDatabaseFile=NULL,enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,
collapseMethod="mean",referenceGeneFile=NULL,referenceGene=NULL,referenceGeneType=NULL,
referenceSet=NULL, minNum=10, maxNum=500,fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,
topThr=10,dNum=20,perNum=1000,lNum=20,is.output=TRUE,outputDirectory=getwd(),
keepGSEAFolder=FALSE,methodType="R",hostName="http://www.webgestalt.org/",
is_parallel=FALSE,nThreads=3)
}

\arguments{
	\item{interestGeneFolder}{
  The folder containing multiple interesting gene files. If \code{enrichMethod} is 
  \code{ORA}, the extension of all files should be \code{txt} and each file can only 
  contain one column: the interesting gene list. If \code{enrichMethod} is \code{GSEA}, 
  the extension of each file should be \code{rnk} and the file should contain 
  two columns separated by tab: the gene list and the corresponding scores.
  }
  \item{interestGeneType}{
	The ID type of the lists in all files. The supported ID type of the WebGestaltR 
	for the selected organism can be found by the function \code{listIDType}. If the 
	\code{organism} is \code{others}, users do not need to set this parameter.
	NOTE: the ID type in all files should be the same.
  }
  \item{enrichMethod}{
	WebGestaltR supports two enrichment analysis methods: 
	ORA (Over-Representation Analysis) and GSEA (Gene Set Enrichment Analysis).
  }
  \item{organism}{
	Currently, WebGestaltR supports 12 organisms. 
	Users can use the function \code{listOrganism} to check the available organisms. 
	Users can also input \code{others} to perform the enrichment analysis based on other 
	organisms not supported by WebGestaltR. For the other organisms, users need to upload 
	the enrichment categories, interesting list and reference list (for ORA method). 
	Because WebGestaltR does not perform the ID mapping for the other organisms, 
	the above uploaded data should have the same ID type.
  }
  \item{enrichDatabase}{
	The functional categories for the enrichment analysis. 
	Users can use the function \code{listGeneset} to check the available functional 
	databases for the selected organism. Users can also input \code{others} to 
	upload the functional database not supported by WebGestaltR for the selected organism. 
  }
  \item{enrichDatabaseFile}{
	If users set \code{organism} as \code{others} or set \code{enrichDatabase} 
	as \code{others}, users need to upload a GMT file as the functional categories 
	for the enrichment analysis. The extension of the file should be \code{gmt} and 
	the first column of the file is the category ID, the second one is the external link 
	for the category. Genes annotated to the category are from the third column. 
	All columns are separated by tab.
  }
  \item{enrichDatabaseType}{
  If users set \code{enrichDatabase} as \code{others}, WebGestaltR will also perform 
  ID mapping for the uploaded GMT file. Thus, users need to set the ID type of the 
  genes in the \code{enrichDatabaseFile}. If users set \code{organism} as \code{others}, 
  users do not need to set this ID type because WebGestaltR will not perform ID mapping 
  for other organisms. The supported ID type of the WebGestaltR for the selected organism 
  can be found by the function \code{listIDType}.
  }
  \item{enrichDatabaseDescriptionFile}{
	Users can also upload a description file for the uploaded \code{enrichDatabaseFile}. 
	The extension of the description file should be \code{des}. The description file 
	contains two columns: the first column is the category ID that should be exactly 
	the same as the category ID in the uploaded \code{enrichDatabaseFile} and the second 
	column is the description of the category. All columns are separated by tab.
  }
  \item{collapseMethod}{
	The method to collapse the duplicate ids for the GSEA method. \code{mean}, 
	\code{median}, \code{min} and \code{max} represent the mean, median, minimum 
	and maximum of scores for the duplicate ids.
  }
  \item{referenceGeneFile}{
	For ORA method, the users need to upload the reference gene list. The extension 
	of the \code{referenceGeneFile} should be \code{txt} and the file can only 
	contain one column: the reference gene list.
  }
  \item{referenceGene}{
	For ORA method, users can also use the R object as the reference gene list.  
	\code{referenceGene} should be an R \code{vector} object containing the 
	reference gene list.
  }
  \item{referenceGeneType}{
	The ID type of the reference gene list. The supported ID type of the 
	WebGestaltR for the selected organism can be found by the function \code{listIDType}. 
	If the \code{organism} is \code{others}, users do not need to set this parameter.
  }
  \item{referenceSet}{
  	Users can directly select the reference set from the existing platform in 
  	the WebGestaltR and do not need to upload the reference set. All existing platform 
  	supported in the WebGestaltR can be found by the function \code{listReferenceSet}. 
  	If \code{referenceGeneFile} and \code{refereneceGene} are \ code{NULL}, WebGestaltR 
  	will use the \code{referenceSet} as the reference gene set. Otherwise, 
  	WebGestaltR will use the user uploaded reference set for the enrichment analysis.
  }
  \item{minNum}{
	WebGestaltR will exclude the categories with the number of annotated genes 
	less than \code{minNum} for the enrichment analysis. The default is \code{10}.
  }
  \item{maxNum}{
	WebGestaltR will exclude the categories with the number of annotated genes 
	larger than \code{maxNum} for the enrichment analysis. The default is \code{500}.
  }
  \item{fdrMethod}{
	For the ORA method, WebGestaltR supports five FDR methods: \code{holm}, 
	\code{hochberg}, \code{hommel}, \code{bonferroni}, \code{BH} and \code{BY}. 
	The default is \code{BH}.
  }
  \item{sigMethod}{
	Two significant methods are available in the WebGestaltR: \code{fdr} and 
	\code{top}. \code{fdr} means the enriched categories are identified based 
	on the FDR and \code{top} means all categories are ranked based on FDR and 
	then selected top categories as the enriched categories. The default is \code{fdr}.
  }
  \item{fdrThr}{
  	The significant level for the \code{fdr} method. The default is \code{0.05}.
  }
  \item{topThr}{
	The threshold for the \code{top} method. The default is \code{10}.
  }
   \item{dNum}{
	The number of enriched categories visualized in the final report. 
	The default is \code{20} and the maximum is \code{100}. A larger \code{dNum} 
	will increase the running time.
  }
   \item{perNum}{
	The number of permutations for the GSEA method. The default is \code{1000}.
  } 
  \item{lNum}{
	The number of categories with the output leading edge genes for the GSEA method. 
	The default is \code{20}.
	\code{Note}: GSEA first ranks the categories based on NES (normalized
	enrichment score) instead of FDR and then outputs the leading edge genes for
	top \code{lNum} categories. Because NES does not necessarily decrease with the 
	increase of the FDR, using \code{sigMethod} defined in WebGestaltR to identify
	the significant categories may cause some categories with outputted leading edge
	genes are not included in the final result even if the number of significant categories 
	is larger than \code{lNum}.
  } 
  \item{is.output}{
  	If \code{is.output} is TRUE, WebGestaltR will create a folder named by the 
  	\code{projectName} and save the mapping results, GO slim summary, enrichment 
  	results and an user-friendly HTML report in the folder. Otherwise, WebGestaltR will 
  	only return an R \code{data.frame} object containing the enrichment results. 
  	If hundreds of gene list need to be analyzed simultaneous, it is better to set 
  	\code{is.output} as FALSE.
  }
  \item{outputDirectory}{
	The output directory for the results.
  } 
  \item{keepGSEAFolder}{
	If \code{keepGSEAFolder} is TRUE, WebGestaltR will keep all folders generated 
	from GSEA tool that contain all figures and tables related to the GSEA analysis.
  }  
  \item{methodType}{
  	For the large ID mapping table (e.g. dbSNP mapping table),
	Users can use \code{R} or \code{Python} function to read it.
	Sometimes \code{Python} code is faster than the \code{R}
	code.
	If users use \code{Python} code to read the mapping
	table, users should first install python and the module \code{pandas}
	in the computer.
  }
  \item{hostName}{
  The server URL for accessing the data. User can use \code{listArchiveURL} function 
  to get all archive version URL.
  }
  \item{is_parallel}{
  If \code{is_parallel} is TRUE, WebGestaltR_batch will use parallel computing to
  simultaneously analyze the lists in all files.
  }
  \item{nThreads}{
  The number of cores used for parallel computing.
  }
}

\value{
	If \code{is.output} is TRUE, each enriched result will be saved in a folder with 
	the name containing the input file name under the \code{outputDirectory}. Otherwise,
	the WebGestaltR_batch function will return a list object containing all results.
	
	If there are errors during the calculation, error message can also be found in the 
	returned list object.
}

\author{
	Jing Wang
}

\examples{
###interestGeneFolder contains multiple .txt files for ORA analysis
refFile<-system.file("extdata","referenceGenes.txt",package="WebGestaltR")
outputDirectory<-getwd()
#enrichResult<-WebGestaltR_batch(interestGeneFolder=interestGeneFolder,
#interestGeneType="genesymbol,enrichMethod="ORA",organism="hsapiens",
#enrichDatabase="pathway_KEGG",referenceGeneFile=refFile, 
#referenceGeneType="genesymbol",is.output=TRUE,
#outputDirectory=outputDirectory,is_parallel=FALSE)

}

\keyword{ methods }
