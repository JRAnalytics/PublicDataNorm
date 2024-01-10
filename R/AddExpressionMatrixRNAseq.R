#' AddExpressionMatrixRNAseq to Meta object
#'
#' @param Metadata Meta object
#' @param Raw.file.path dir path in which the raw single cell exrpession is saved
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
AddExpressionMatrixRNAseq <- function(Metadata=NULL,
                                      Raw.file.path,
                                      name,
                                      Export=T,
                                      name.local.file = NULL,
                                      force.replace=F ) {

   Omics.type = "RNAseq"

   path = Raw.file.path

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


      message("Local import")

      l <-length(names(Metadata))
      lf <- list.files(path)

      if(length(lf)>1){print(c(message("There is more than one files in Dir :"),lf))}

      if(all(str_detect(lf, ".rds|.txt|.csv|.tsv|.mtx", negate = FALSE)==F)){stop("No '*.rds' or '*.txt' or '*.csv' or '*.mtx' files in set directory. \n change path or add file")}

      if(!is.null(name.local.file)) {

        filepath <- paste(path,name.local.file,sep="/")
        message(paste("Loading", name.local.file, "file"))
        if(str_detect(name.local.file, ".rds", negate = FALSE)){
          dt <- readRDS(filepath)}
        else {
          if(str_detect(name.local.file, ".txt", negate = FALSE)){

            dt <- suppressWarnings(as.data.frame(data.table::fread(filepath)))

            if(length(dt[,1])!=length(unique( dt[,1]))) { } else {  rownames(dt) <- dt[,1]}


          }  else {
            if(str_detect(name.local.file, ".csv", negate = FALSE)){
              dt <- suppressWarnings(as.data.frame(data.table::fread(filepath)))
              if(length(dt[,1])!=length(unique( dt[,1]))) { } else {  rownames(dt) <- dt[,1]}
            }else {

              if(str_detect(name.local.file, ".tsv", negate = FALSE)){
                dt <- suppressWarnings(as.data.frame(data.table::fread(filepath)))
                if(length(dt[,1])!=length(unique( dt[,1]))) { } else {  rownames(dt) <- dt[,1]}
              } else {


              }# tsv


            }#csv
          }#txt
        }#rds




        if(length(Metadata)>=1) {

          if(!all(str_detect(names(Metadata),paste0(name,".matrix")))==F){
            message("An Object with the same name already exist in MetaObject")
            if(force.replace==F){stop("set force.replace==T to subset object.")}
            message("Subsetting object.")
            Metadata[[paste0(name,".matrix")]] <- dt    } else { Metadata[[l+1]] <- dt
            names(Metadata)[l+1] <- paste0(name,".matrix")}

          if(length(attributes(Metadata)$Data.Type)<length(Metadata)){
            attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type, "Count")

            if(Export==T){attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes") } else {attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"No") }
          }

          return(Metadata)} # Metadat >1
        else {

          Metadata$mat <-  dt
          names(Metadata)[1] <- paste0(name,".matrix")

          if(length(attributes(Metadata)$Data.Type)<length(Metadata)){
            attributes(Metadata)$Data.Type <-  c("Count")
            if(Export==T){attributes(Metadata)$Export <- "Yes" } else {attributes(Metadata)$Export <- "No" }}
          return(Metadata)
        } # Metadata = 0



      } else {


        if(length(lf[str_detect(lf, "matrix")])>1) {

          for (i in lf[str_detect(lf, "matrix")]) {
            l <- length(Metadata)
            message(paste("Loading", i, "file"))

            if(str_detect(i, ".rds", negate = FALSE)){Metadata[[1]] <- readRDS(i)} else {
              if(str_detect(i, ".txt", negate = FALSE)){

                dt <- suppressWarnings(as.data.frame(data.table::fread(i)))
                rownames(dt) <- dt[,1]
                if(!all(str_detect(names(Metadata),paste0(name,".matrix")))==F){
                  message("An Object with the same name already exist in MetaObject")
                  if(force.replace==F){stop("set force.replace==T to subset object.")}
                  message("Subsetting object.")
                  Metadata[[paste0(name,".matrix")]] <- dt    } else { Metadata[[l+1]] <- dt
                  names(Metadata)[l+1] <- paste0(name,".matrix")}

                attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type, "Count")
                if(Export==T){attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes") } else {attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"No") }
              }  else {
                if(str_detect(i, ".csv", negate = FALSE)){
                  dt <- suppressWarnings(as.data.frame(data.table::fread(i)))
                  rownames(dt) <- dt[,1]
                  if(!all(str_detect(names(Metadata),paste0(name,".matrix")))==F){
                    message("An Object with the same name already exist in MetaObject")
                    if(force.replace==F){stop("set force.replace==T to subset object.")}
                    message("Subsetting object.")
                    Metadata[[paste0(name,".matrix")]] <- dt    } else { Metadata[[l+1]] <- dt
                    names(Metadata)[l+1] <- paste0(name,".matrix")}

                  attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type, "Count")
                  if(Export==T){attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes") } else {attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"No") }} else {

                    if(str_detect(i, ".tsv", negate = FALSE)){
                      dt <- suppressWarnings(as.data.frame(data.table::fread(i)))
                      rownames(dt) <- dt[,1]
                      if(!all(str_detect(names(Metadata),paste0(name,".matrix")))==F){
                        message("An Object with the same name already exist in MetaObject")
                        if(force.replace==F){stop("set force.replace==T to subset object.")}
                        message("Subsetting object.")
                        Metadata[[paste0(name,".matrix")]] <- dt    } else { Metadata[[l+1]] <- dt
                        names(Metadata)[l+1] <- paste0(name,".matrix")}

                    }
                    attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type, "Count")
                    if(Export==T){attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes") } else {attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"No") }}
              } }

            names(Metadata)[l+1] <- paste0(name,".matrix.",which(lf[str_detect(lf, "matrix")]%in%i))

          } } else {

            if(length(Metadata)>=1) {


              lf <- lf[str_detect(lf, "matrix")]
              message(paste("Loading", lf, "file"))

              if(str_detect(lf, ".rds", negate = FALSE)){Metadata[[l+1]] <- readRDS(lf)} else {

                if(str_detect(lf, ".txt", negate = FALSE)){
                  dt <- suppressWarnings(as.data.frame(data.table::fread(lf)))
                  rownames(dt) <- dt[,1]
                  dt <- dt[,colnames(Metadata[[1]])]

                  if(!all(str_detect(names(Metadata),paste0(name,".matrix")))==F){
                    message("An Object with the same name already exist in MetaObject")
                    if(force.replace==F){stop("set force.replace==T to subset object.")}
                    message("Subsetting object.")
                    Metadata[[paste0(name,".matrix")]] <- dt    } else { Metadata[[l+1]] <- dt
                    names(Metadata)[l+1] <- paste0(name,".matrix")}

                  if(length(attributes(Metadata)$Data.Type)<length(Metadata)){
                    attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type, "Count")
                    if(Export==T){attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes") } else {attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"No") }}  else {
                    }
                  if(str_detect(lf, ".csv", negate = FALSE)){
                    dt <- suppressWarnings(as.data.frame(data.table::fread(lf)))
                    rownames(dt) <- dt[,1]
                    dt <- dt[,colnames(Metadata[[1]])]

                    if(!all(str_detect(names(Metadata),paste0(name,".matrix")))==F){
                      message("An Object with the same name already exist in MetaObject")
                      if(force.replace==F){stop("set force.replace==T to subset object.")}
                      message("Subsetting object.")
                      Metadata[[paste0(name,".matrix")]] <- dt    } else { Metadata[[l+1]] <- dt
                      names(Metadata)[l+1] <- paste0(name,".matrix")}

                    if(length(attributes(Metadata)$Data.Type)<length(Metadata)){
                      attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type, "Count")
                      if(Export==T){attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes") } else {attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"No") }
                    }} else {

                      if(str_detect(lf, ".tsv", negate = FALSE)){
                        dt <- suppressWarnings(as.data.frame(data.table::fread(lf)))
                        rownames(dt) <- dt[,1]
                        dt <- dt[,colnames(Metadata[[1]])]

                        if(!all(str_detect(names(Metadata),paste0(name,".matrix")))==F){
                          message("An Object with the same name already exist in MetaObject")
                          if(force.replace==F){stop("set force.replace==T to subset object.")}
                          message("Subsetting object.")
                          Metadata[[paste0(name,".matrix")]] <- dt    } else { Metadata[[l+1]] <- dt
                          names(Metadata)[l+1] <- paste0(name,".matrix")}


                        if(length(attributes(Metadata)$Data.Type)<length(Metadata)){
                          attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type, "Count")
                          if(Export==T){attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes") } else {attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"No") }
                        }
                      }
                    }
                }

              }

            }

            else {

              lf <- lf[str_detect(lf, "matrix")]
              message(paste("Loading", lf, "file"))
              if(str_detect(lf, ".rds", negate = FALSE)){Metadata[[1]] <- readRDS(lf)} else {
                if(str_detect(lf, ".txt", negate = FALSE)){

                  dt <- suppressWarnings(as.data.frame(data.table::fread(lf)))
                  rownames(dt) <- dt[,1]

                  Metadata[[1]] <- dt
                  attributes(Metadata)$Data.Type <-  c("Count")
                  if(Export==T){attributes(Metadata)$Export <- "Yes" } else {attributes(Metadata)$Export <- "No" } }  else {
                    if(str_detect(lf, ".csv", negate = FALSE)){
                      dt <- suppressWarnings(as.data.frame(data.table::fread(lf)))
                      rownames(dt) <- dt[,1]

                      Metadata[[1]] <- dt
                      attributes(Metadata)$Data.Type <-  c("Count")
                      if(Export==T){attributes(Metadata)$Export <- "Yes" } else {attributes(Metadata)$Export <- "No" } }else {

                        if(str_detect(lf, ".tsv", negate = FALSE)){
                          dt <- suppressWarnings(as.data.frame(data.table::fread(lf)))
                          rownames(dt) <- dt[,1]

                          Metadata[[1]] <- dt
                          attributes(Metadata)$Data.Type <-  c("Count")
                          if(Export==T){attributes(Metadata)$Export <- "Yes" } else {attributes(Metadata)$Export <- "No" } }}
                  }
              }

              names(Metadata)[1] <- paste0(name,".matrix")




            }}


        ###marche pas le readme!


        return(Metadata)


      }

  }else{ stop("No meta data object found")} #else meta data object not found











} # function
