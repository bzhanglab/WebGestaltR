\name{WebGestaltR}
\alias{WebGestaltR}

\title{
	Comprehensive R function for the enrichment analysis
}

\description{
	WebGestaltR function can perform two popular enrichment analyses: 
	ORA (Over-Representation Analysis) and GSEA (Gene Set Enrichment Analysis). 
	Based on the user uploaded gene list or gene list with scores (for GSEA method),  
	WebGestaltR function will first map the gene list to the entrez gene ids and then 
	summary the gene list based on the GO (Gene Ontology) Slim. After performing the 
	enrichment analysis, WebGestaltR function also returns an user-friendly HTML report 
	containing the ID mapping table, GO Slim summary result and the enrichment analysis 
	result. If the functional categories have the DAG (directed acyclic graph) structure, 
	the structure of the enriched categories can also be visualized in the report.
}

\usage{
WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
enrichDatabase="geneontology_Biological_Process",enrichDatabaseFile=NULL, 
enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL, 
interestGene=NULL,interestGeneType=NULL,collapseMethod="mean",referenceGeneFile=NULL,
referenceGene=NULL,referenceGeneType=NULL,referenceSet=NULL, minNum=10, maxNum=500,
fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,topThr=10,dNum=20,perNum=1000,
lNum=20,is.output=TRUE,outputDirectory=getwd(),projectName=NULL,keepGSEAFolder=FALSE,
methodType="R",dagColor="continuous",hostName="http://www.webgestalt.org/")
}

\arguments{
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
  \item{interestGeneFile}{
  If \code{enrichMethod} is \code{ORA}, the extension of the \code{interestGeneFile} 
  should be \code{txt} and the file can only contain one column: the interesting gene 
  list. If \code{enrichMethod} is \code{GSEA}, the extension of the 
  \code{interestGeneFile} should be \code{rnk} and the file should contain 
  two columns separated by tab: the gene list and the corresponding scores. 
  }
  \item{interestGene}{
	Users can also use the R object as the input. If \code{enrichMethod} is \code{ORA}, 
	\code{interestGene} should be an R \code{vector} object containing the interesting 
	gene list. If \code{enrichMethod} is \code{GSEA}, \code{interestGene} should be an R 
	\code{data.frame} object containing two columns: the gene list and the corresponding 
	scores.
  }
  \item{interestGeneType}{
	The ID type of the interesting gene list. The supported ID type of the WebGestaltR 
	for the selected organism can be found by the function \code{listIDType}. If the 
	\code{organism} is \code{others}, users do not need to set this parameter.
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
  \item{projectName}{
  	The name of the project. If \code{projectName} is Null, WebGestaltR will 
  	use time stamp as the project name.
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
   \item{dagColor}{
   If \code{dagColor} is binary, the significant terms in the DAG structure will be colored by
   red for ORA method or red (positive related) and blue (negative related) for GSEA
   method. If \code{dagColor} is continous, the significant terms in the DAG structure will
   be colored by the red gradient for ORA method or red (positive related) and blue 
   (negative related) gradient for GSEA method.based on the corresponding FDR.
  }
  \item{hostName}{
  The server URL for accessing the data. User can use \code{listArchiveURL} function 
  to get all archive version URL.
  }
}

\value{
	The WebGestaltR function not only outputs the user-friendly HTML report containing 
	the ID mapping table, GO Slim summary result and the enrichment analysis result but 
	also outputs an R object containing the enrichment analysis result.
}

\note{
	Because WebGestaltR will read the data from the server for the ID mapping and GO 
	Slim summary, the running time for the WebGestaltR function will be also based on 
	the internet speed. Generally, it will take around one minute to perform the whole 
	analysis. Because of the huge number of the Gene ontology categories, running GSEA 
	analysis for these categories may take one to six minutes based on the different 
	size of the uploaded gene rank list. Decreasing the parameter \code{maxNum} can 
	reduce the running time.
	
	If the user has any problem for the Network Visualization in the HTML report, please
	follow the instruction in the  http://cytoscapeweb.cytoscape.org/tutorial.
}

\author{
	Jing Wang
}

\examples{

#######ORA example#########
geneFile<-system.file("extdata","interestingGenes.txt",package="WebGestaltR")
refFile<-system.file("extdata","referenceGenes.txt",package="WebGestaltR")
outputDirectory<-getwd()
#enrichResult<-WebGestaltR(enrichMethod="ORA",organism="hsapiens",
#enrichDatabase="pathway_KEGG",interestGeneFile=geneFile, 
#interestGeneType="genesymbol",referenceGeneFile=refFile, 
#referenceGeneType="genesymbol",is.output=TRUE,
#outputDirectory=outputDirectory,projectName=NULL)
	
#######GSEA example#########
#rankFile<-system.file("extdata","GeneRankList.rnk",package="WebGestaltR")
#outputDirectory<-getwd()
#enrichResult<-WebGestaltR(enrichMethod="GSEA",organism="hsapiens",
#enrichDatabase="pathway_KEGG",interestGeneFile=rankFile,
#interestGeneType="genesymbol", collapseMethod="mean",
#is.output=TRUE,outputDirectory=outputDirectory)
}

\keyword{ methods }
