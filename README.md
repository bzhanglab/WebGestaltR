# WebGestalt

WebGestalt R package is the R version of our well-known web application tool WebGestalt (www.webgestalt.org) that has on average 27,000 users from 140 countries and territories per year and has also been cited 371 in 2016. The advantage of this R package is that it can be easily integrated to other pipelines or simultaneously analyze multiple gene lists.

WebGestaltR function can perform popular enrichment analyses: ORA (Over-Representation Analysis), GSEA (Gene Set Enrichment Analysis) and NTA (Network Topology Analysis). Based on the user-uploaded gene list or gene list with scores (for GSEA method), WebGestaltR function will first map the gene list to entrez gene IDs and then summarize the gene list based on the GO (Gene Ontology) Slim data. After performing the enrichment analysis, WebGestaltR function also returns an user-friendly HTML report containing GO Slim summary and enrichment analysis result. If the functional categories have the DAG (directed acyclic graph) structure, the structure of the enriched categories can also be visualized in the report.

## Installation

Since WebGestaltR v2.0.0, Rust is used for core computations in the R package. Therefore, to install WebGestaltR, please download and install Rust from [https://www.rust-lang.org/learn/get-started](https://www.rust-lang.org/learn/get-started). For Mac, Linux, or Unix users, Rust can be installed from the command line, and Windows users can download a GUI installer.
