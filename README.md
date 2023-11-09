# WebGestaltR

[![R-CMD-check](https://github.com/iblacksand/WebGestaltR/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/iblacksand/WebGestaltR/actions/workflows/check-standard.yaml)

> [!IMPORTANT]
> The new version of WebGesaltR requires Rust, which must be installed on your device prior to installing or updating the package from CRAN. See the installation section for more information.

WebGestalt R package is the R version of our well-known web application tool WebGestalt (www.webgestalt.org) that has on average 27,000 users from 140 countries and territories per year and has also been cited 371 in 2016. The advantage of this R package is that it can be easily integrated to other pipelines or simultaneously analyze multiple gene lists.

WebGestaltR function can perform popular enrichment analyses: ORA (Over-Representation Analysis), GSEA (Gene Set Enrichment Analysis) and NTA (Network Topology Analysis). Based on the user-uploaded gene list or gene list with scores (for GSEA method), WebGestaltR function will first map the gene list to entrez gene IDs and then summarize the gene list based on the GO (Gene Ontology) Slim data. After performing the enrichment analysis, WebGestaltR function also returns an user-friendly HTML report containing GO Slim summary and enrichment analysis result. If the functional categories have the DAG (directed acyclic graph) structure, the structure of the enriched categories can also be visualized in the report.

## Installation

---

### Requirements

- R (>= 4.0.0)
- Rust (>= 1.63.0)

---

Since WebGestaltR v1.0.0, Rust is used for core computations in the R package. Therefore, to install WebGestaltR, please download and install Rust from [https://www.rust-lang.org/tools/install](https://www.rust-lang.org/tools/install). For Mac, Linux, or Unix users, Rust can be installed from the command line, and Windows users can download a GUI installer.

Make sure you restart your terminal after installing Rust to ensure the Rust compiler is available in your path. You can check that Rust is installed correctly by running `rustc --version` in your terminal.

After installing Rust, you can install WebGestaltR by running the following command in an R session:

```R
# install.packages("devtools") # run if devtools not already installed
devtools::install_github("iblacksand/WebGestaltR")
```

During installation, the Rust compiler will be called to build the computation library used by WebGestaltR. If you run into problems with installation of the new version, please [open a new issue](https://github.com/iblacksand/WebGestaltR/issues/new/choose).

## Changes

> [!NOTE]
> Besides the change in installation, there should be no difference in how the R package performs for existing use-cases. If you experience any difference in results that are not due to the data-update, that is considered a bug. [Please report the changes you experience in a new issue](https://github.com/iblacksand/WebGestaltR/issues/new/choose).

WebGestaltR's core was re-written in Rust, which dramatically increased performance, with up to 15x the speed of previous versions. The new version also supports metabolomics, with support for 15 different ID types.
