% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TCGAbuild.R
\name{TCGA.build}
\alias{TCGA.build}
\title{TCGA.build allow to build a gene expression data frame imported from the TCGA sites}
\usage{
TCGA.build(
  ask,
  files,
  cases,
  experimental.strategy,
  merging.col = c("gene_id"),
  legacy = F
)
}
\arguments{
\item{ask}{a character "Raw"  = raw count, "TPM" = tpm normalisation gene expression, "FKPM" = fkpm normalisation gene expression if it exist}

\item{files}{a character list of files download from the TCGA in sub directory}

\item{cases}{sample id}

\item{experimental.strategy}{a character "Gene expression array" or "RNA-Seq"}

\item{merging.col}{column name to merge files c("gene_id", "gene_name")}

\item{legacy}{True or False. True : Leagcy files in TCGA}
}
\value{
a data frame with gene in in rows as Ensembl, and samples count or normalized count in column
}
\description{
TCGA.build allow to build a gene expression data frame imported from the TCGA sites
}
\examples{

"none"

}
