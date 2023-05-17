#' TCGA.build allow to build a gene expression data frame imported from the TCGA sites
#'
#' @param ask a character "Raw"  = raw count, "TPM" = tpm normalisation gene expression, "FKPM" = fkpm normalisation gene expression if it exist
#' @param files a character list of files download from the TCGA in sub directory
#' @param cases sample id
#' @param genome version of the human genome "hg19" or "hg38"
#' @param experimental.strategy a character "Gene expression array" or "RNA-Seq"
#' @param merging.col column name to merge files c("gene_id", "gene_name")
#' @param legacy True or False. True : Leagcy files in TCGA
#' @import data.table
#' @import TCGAbiolinks
#' @import dplyr
#' @import DT
#' @importFrom plyr alply
#' @importFrom purrr reduce
#' @return a data frame with gene in in rows as Ensembl, and samples count or normalized count in column
#' @examples
#'
#' "none"
#'
#' @export
#'
#'
#'
TCGA.build <- function(ask,
  files,
  cases,
  genome = "hg38",
  experimental.strategy,
  merging.col = c("gene_id"),
  legacy = F
){
  skip <- unique((ifelse(experimental.strategy == "Gene expression array",1,0)))

  if(length(skip) > 1) stop("It is not possible to handle those different platforms together")


  ret <- plyr::alply(
    .data = seq_along(files),
    .margins = 1,
    .fun = function(i,cases){

      data <- as.data.frame(fread(
        input = files[i],
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        skip = skip
      ))


  if(legacy == F){
      if(ask=="Raw"){ data <- data[-c(1:4),c(merging.col,"unstranded")]}
      if(ask=="TPM"){ data <- data[-c(1:4),c(merging.col,"tpm_unstranded")]}
      if(ask=="FKPM"){ data <- data[-c(1:4),c(merging.col,"fpkm_unstranded")]}
    if(ask=="FKPM_UQ"){ data <- data[-c(1:4),c(merging.col,"fpkm_uq_unstranded")]}}

  if(legacy == T){

    if(ask=="Raw"){ data <- data[,c(merging.col,"raw_count")]}
    if(ask=="Normalized results"){ data <- data[,c(merging.col,"normalized_count")]}

  }
      if(!missing(cases)) {
        #assay.list <<- gsub(" |\\(|\\)|\\/","_",colnames(data)[2:ncol(data)])
        # We will use this because there might be more than one col for each samples
        setnames(data,colnames(data)[2:ncol(data)],
                 paste0(gsub(" |\\(|\\)|\\/","_",colnames(data)[2:ncol(data)]),"_",cases[i]))
      }
      data
    },.progress = "time",cases = cases)



  df <- purrr::reduce(
    ret,
    dplyr::full_join,
    by = merging.col
  )
  df <- as.data.frame(df)
  rownames(df) <- df[,merging.col]
  df[,merging.col] <- NULL



  patient <- query$results[[1]]$cases.submitter_id
  if(length(patient)> length(unique(patient))){   patient <- query$results[[1]]$sample.submitter_id }
  colnames(df) <- patient


  return(df)
}
