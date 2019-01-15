### R code from vignette source 'WebGestaltR.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Example
###################################################

library("WebGestaltR")
#######ORA example#########
#interestGeneFile <- system.file("extdata","interestingGenes.txt",package="WebGestaltR")
#referenceGeneFile <- system.file("extdata","referenceGenes.txt",package="WebGestaltR")
#outputDirectory <- getwd()
#enrichResult <-WebGestaltR(enrichMethod="ORA", organism="hsapiens", enrichDatabase="geneontology_Biological_Process", enrichDatabaseFile=NULL, enrichDatabaseType=NULL, enrichDatabaseDescriptionFile=NULL, interestGeneFile=interestGeneFile, interestGene=NULL, interestGeneType="genesymbol", collapseMethod="mean", referenceGeneFile=referenceGeneFile, referenceGene=NULL, referenceGeneType="genesymbol", referenceSet=NULL, minNum=10, maxNum=500, fdrMethod="BH", sigMethod="fdr", fdrThr=0.01, topThr=10, reportNum=20, perNum=1000, lNum=20, is.output=TRUE, outputDirectory=outputDirectory, projectName=NULL, keepGSEAFolder=FALSE,methodType="R",dagColor="continuous",hostName="http://www.webgestalt.org/")

########GSEA example#######
#geneRankFile <- system.file("extdata","GeneRankList.rnk",package="WebGestaltR")
#outputDirectory <- getwd()
#enrichResult <-WebGestaltR(enrichMethod="GSEA", organism="hsapiens", enrichDatabase="pathway_KEGG", enrichDatabaseFile=NULL, enrichDatabaseType=NULL, enrichDatabaseDescriptionFile=NULL, interestGeneFile=geneRankFile, interestGene=NULL, interestGeneType="genesymbol", collapseMethod="mean", referenceGeneFile=NULL, referenceGene=NULL, referenceGeneType=NULL, referenceSet=NULL, minNum=10, maxNum=500, fdrMethod="BH", sigMethod="fdr", fdrThr=0.01, topThr=10, reportNum=20, perNum=1000, lNum=20, is.output=TRUE, outputDirectory=outputDirectory, projectName=NULL, keepGSEAFolder=FALSE,methodType="R",hostName="http://www.webgestalt.org/")



###################################################
### code chunk number 2: Example
###################################################
#library("WebGestaltR")
#url <- listArchiveURL()


###################################################
### code chunk number 3: Example
###################################################
#library("WebGestaltR")
#geneFile<-system.file("extdata","interestingGenes.txt",package="WebGestaltR")
#interestGene <- formatCheck(dataType="list",inputGeneFile=geneFile,inputGene=NULL)


###################################################
### code chunk number 4: Example
###################################################
#library("WebGestaltR")
#interestGeneFile <- system.file("extdata","interestingGenes.txt",package="WebGestaltR")
#idmap <- IDMapping(organism="hsapiens", dataType="list", inputGeneFile=interestGeneFile, inputGene=NULL, sourceIdType, targetIdType, collapseMethod="mean", mappingOutput=FALSE, outputFileName="", methodType="R", hostName="http://www.webgestalt.org/")


###################################################
### code chunk number 5: Example
###################################################
#library("WebGestaltR")
#geneListFile <- system.file("extdata","GOSlimExample.txt",package="WebGestaltR")
#geneList <- read.table(geneListFile, header=FALSE, sep="\t", stringsAsFactors=FALSE)
#geneList <- as.vector(as.matrix(geneList))
#outputFile <- paste(getwd(),"/GOSlimSummary",sep="")
#GOSlimSummary(organism="hsapiens", genelist=geneList, outputFile=outputFile, outputType="pdf", hostName="http://www.webgestalt.org/")


###################################################
### code chunk number 6: Example
###################################################
#library("WebGestaltR")
#organism <- listOrganism(hostName="http://www.webgestalt.org/")


###################################################
### code chunk number 7: Example
###################################################
#library("WebGestaltR")
#geneSet <- listGeneSet(organism="hsapiens",hostName="http://www.webgestalt.org/")


###################################################
### code chunk number 8: Example
###################################################
#ibrary("WebGestaltR")
#idType <- listIDType(organism="hsapiens", hostName="http://www.webgestalt.org/")


###################################################
### code chunk number 9: Example
###################################################
#library("WebGestaltR")
#referenceSet <- listReferenceSet(organism="hsapiens", hostName="http://www.webgestalt.org/")


