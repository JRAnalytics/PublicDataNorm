#' AddExpressionMatrixTCGA from TCGA to Meta object
#'
#' @param Metadata Meta object
#' @param query query file from TCGAimportExpression function
#' @param data.norm a character "Raw" ,"TPM", "FKPM"
#' @param Raw.file.path dir path in which the GDC project is saved, or local files are saved
#' @param name if loca=True, names to apply in Metadata object slot
#' @param name.local.file name file of interest in path directory
#' @param force.replace set as F. T : replace an already object with the same name
#' @param Export  TRUE or FALSE. If data to be Exported, set T.
#' @importFrom utils menu
#' @importFrom Matrix readMM
#' @import data.table
#' @return a data.frame in the Meta Object
#' @export
#'
#' @examples "none"
AddExpressionMatrixTCGA <- function(Metadata=NULL,
                                    query,
                                    data.norm,
                                    Raw.file.path,
                                    name,
                                    Export=T,
                                    name.local.file = NULL,
                                    force.replace=F ) {

  Omics.type = "RNAseq"
  path =Raw.file.path


  if(is.null(Metadata)){
    Metadata = list()
    attributes(Metadata)$Omics.type = Omics.type
  }

  if(!is.null(Metadata)){

    if(!attributes(Metadata)$Omics.type==Omics.type){
      warning(paste("Omics.type is", Omics.type,"different than attributes(Metadata)$Omics.type",
                    attributes(Metadata)$Omics.type, "\nWill be replace"))

      attributes(Metadata)$Omics.type = Omics.type

    }

    if(!is.list(Metadata)){
      stop("Metadata should be a list.")}

    setwd(path)
        message("Fetching data from TCGA portal")
        query <- query

        if(is.null(query)){message("Query data form adding expression matrix to Meta object needed \n launch : query <- GDCquery() with associated project")}


        project <- query$project

        source <- "harmonized"
        files <- file.path(
          project, source,
          gsub(" ","_",query$results[[1]]$data_category),
          gsub(" ","_",query$results[[1]]$data_type),
          gsub(" ","_",query$results[[1]]$file_id),
          gsub(" ","_",query$results[[1]]$Rename)
        )


        if(dir.exists(file.path("GDCdata",project, source))){message(paste(project, "data found"))}

        files <- file.path("GDCdata", files)


        cases <- ifelse(grepl("TCGA|TARGET",query$results[[1]]$project %>% unlist()),query$results[[1]]$cases,query$results[[1]]$sample.submitter_id)


        message("Building gene expression data")

        y <- TCGA.build(
          ask = data.norm,
          files = files,
          cases = cases,
          experimental.strategy = unique(query$results[[1]]$experimental_strategy))



        if(!all(colnames(Metadata[[1]])==colnames(y))){stop("Samples are not the sames accross data matrix expression")}

        y <- y[,colnames(Metadata[[1]])]

        l <-length(names(Metadata))

        Metadata[[l+1]]<- y
        names(Metadata)[l+1] <- c(paste0(data.norm,".",project,".matrix"))
        if(length(attributes(Metadata)$Data.Type)<length(Metadata)){
          if(l==0) {   attributes(Metadata)$Data.Type <-  c("Expression.Matrix")
          if(Export==T){attributes(Metadata)$Export <- "Yes" } else {attributes(Metadata)$Export <- "No" }
          }} else {
            if(length(attributes(Metadata)$Data.Type)<length(Metadata)){
              attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"Expression.Matrix")

              if(Export==T){attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes") } else {attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"No") } }
          }



    return(Metadata)





  }else{ stop("No meta data object found")} #else meta data object not found











} # function
