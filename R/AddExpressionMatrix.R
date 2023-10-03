#' AddExpressionMatrix to Meta object
#'
#' @param Metadata Meta object
#' @param query query file from TCGAimportExpression function
#' @param data.norm a character "Raw" ,"TPM", "FKPM"
#' @param path dir path in which the GDC project is saved, or local files are saved
#' @param local T of F. If F, use query form. If F, add expression patrice from local file
#' @param name if loca=True, names to apply in Metadata object slot
#' @param Omics.type one of this category : "RNAseq", "Single.Cell", "Microarray", "Spatial".
#' @param name.local.file name file of interest in path directory
#' @param Cell.file name of single cell cell annotation file specific to expression matrix. Not mandatory.
#' @param Genes.file name of single cell gene annotation file specific to expression matrix. Not mandatory.
#' @param force.replace set as F. T : replace an already object with the same name
#' @param Raw TRUE or FALSE. If Raw data, to be specified.
#' @importFrom utils menu
#' @importFrom Matrix readMM
#' @import data.table
#' @return a data.frame in the Meta Object
#' @export
#'
#' @examples "none"
AddExpressionMatrix <- function(Metadata=NULL, local = c(T, F) , Omics.type = c("RNAseq", "Single.Cell", "Microarray", "Spatial"),
                                Cell.file=NULL, Genes.file=NULL, query, data.norm, path, name, Raw=T, name.local.file = NULL, force.replace=F ) {

  if(is.null(Omics.type)){
    if(!inherits(Omics.type, "character")){ stop("Omics.type  is not a character string")}

    stop("Omics.type must be from one of this character c(RNAseq, Single.Cell, Microarray, Spatial)")}

  if(is.null(Metadata)){
    Metadata = list()
    attributes(Metadata)$Omics.type = Omics.type
    }

  if(!is.null(Metadata)){

    if(!attributes(Metadata)$Omics.type==Omics.type){
      warning(paste("Omics.type is", Omics.type,"different than attributes(Metadata)Omics.type",
                    attributes(Metadata)$Omics.type, "\nWill be replace"))

              attributes(Metadata)$Omics.type = Omics.type

              }

    if(!is.list(Metadata)){
      stop("Metadata should be a list.")}

  setwd(path)
  if(local== T){


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

               if(str_detect(name.local.file, ".mtx", negate = FALSE)){

                 dt <- readMM(filepath)
                 dt = as.matrix(dt)

                 if(is.null(colnames(dt))){
                   message(paste(name.local.file,"has no colnames. A Cell.csv file may be associated in raw data directory."))

                   if(!is.null(Cell.file)){
                     message(paste("Loading",Cell.file ))
                     Cells <- read.csv(file.path(path,Cell.file))
                     if("cell_name"%in%colnames(Cells)){
                     rownames(Cells)  = Cells$cell_name} else { stop("Cell.file mus have a colnames specified 'cell_name'")}


                     if(length(rownames(Cells)==length(colnames(dt)))) {  colnames(dt) = rownames(Cells)} else
                     { stop(paste(Cell.file, "has not the same number of cells than column of expression matrix."))
                         }



                     }
                 }


                 if(is.null(rownames(dt))){
                   message(paste(name.local.file,"has no rownames A Genes.csv file may be associated in raw data directory."))

                   if(!is.null(Genes.file)){
                     message(paste("Loading",Genes.file ))
                     Genes <- read.table(file.path(path,Genes.file), quote="\"", comment.char="")
                     if(length(Genes[,1])==nrow(dt)){
                    Genes  = Genes[,1]

                     rownames(dt) =  Genes
                     } else { stop(paste(Genes.file, "has not the same length as rows of expression matrix"))}}
                 }

               } # mtx
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
        attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type, "Expression.Matrix")

        if(Raw==T){attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"Yes") } else {attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"No") }
}

        return(Metadata)} # Metadat >1
      else {

        Metadata$mat <-  dt
        names(Metadata)[1] <- paste0(name,".matrix")

        if(length(attributes(Metadata)$Data.Type)<length(Metadata)){
        attributes(Metadata)$Data.Type <-  c("Expression.Matrix")
        if(Raw==T){attributes(Metadata)$Raw.data <- "Yes" } else {attributes(Metadata)$Raw.data <- "No" }}
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

            attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type, "Expression.Matrix")
            if(Raw==T){attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"Yes") } else {attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"No") }
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

                attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type, "Expression.Matrix")
                if(Raw==T){attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"Yes") } else {attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"No") }} else {

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
                    attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type, "Expression.Matrix")
                    if(Raw==T){attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"Yes") } else {attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"No") }}
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
        attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type, "Expression.Matrix")
        if(Raw==T){attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"Yes") } else {attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"No") }}  else {
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
            attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type, "Expression.Matrix")
            if(Raw==T){attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"Yes") } else {attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"No") }
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
              attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type, "Expression.Matrix")
              if(Raw==T){attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"Yes") } else {attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"No") }
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
        attributes(Metadata)$Data.Type <-  c("Expression.Matrix")
        if(Raw==T){attributes(Metadata)$Raw.data <- "Yes" } else {attributes(Metadata)$Raw.data <- "No" } }  else {
          if(str_detect(lf, ".csv", negate = FALSE)){
            dt <- suppressWarnings(as.data.frame(data.table::fread(lf)))
            rownames(dt) <- dt[,1]

            Metadata[[1]] <- dt
            attributes(Metadata)$Data.Type <-  c("Expression.Matrix")
            if(Raw==T){attributes(Metadata)$Raw.data <- "Yes" } else {attributes(Metadata)$Raw.data <- "No" } }else {

              if(str_detect(lf, ".tsv", negate = FALSE)){
                dt <- suppressWarnings(as.data.frame(data.table::fread(lf)))
                rownames(dt) <- dt[,1]

                Metadata[[1]] <- dt
                attributes(Metadata)$Data.Type <-  c("Expression.Matrix")
                if(Raw==T){attributes(Metadata)$Raw.data <- "Yes" } else {attributes(Metadata)$Raw.data <- "No" } }}
        }
        }

      names(Metadata)[1] <- paste0(name,".matrix")




      }}


###marche pas le readme!


    return(Metadata)


  }} else { # if local ==T
    message("Fetching data from TCGA portal")
  query <- query

  if(is.null(query)){message("Query data form adding expression matrix to Meta object needed \n launch : query <<- GDCquery() with associated project")}


  project <- query$results[[1]]$project

  source <- ifelse(query$legacy,"legacy","harmonized")
  files <- file.path(
    project, source,
    gsub(" ","_",query$results[[1]]$data_category),
    gsub(" ","_",query$results[[1]]$data_type),
    gsub(" ","_",query$results[[1]]$file_id),
    gsub(" ","_",query$results[[1]]$file_name)
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

  y <- y[,colnames(Metadata[[1]])]

  if(!all(colnames(Metadata[[1]])==colnames(y))){stop("Samples are not the sames accross data matrix expression")}


  l <-length(names(Metadata))

  Metadata[[l+1]]<- y
  names(Metadata)[l+1] <- c(paste0(data.norm,".",project,".matrix"))
  if(length(attributes(Metadata)$Data.Type)<length(Metadata)){
  if(l==0) {   attributes(Metadata)$Data.Type <-  c("Expression.Matrix")
  if(Raw==T){attributes(Metadata)$Raw.data <- "Yes" } else {attributes(Metadata)$Raw.data <- "No" }
  }} else {
    if(length(attributes(Metadata)$Data.Type)<length(Metadata)){
      attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"Expression.Matrix")

  if(Raw==T){attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"Yes") } else {attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"No") } }
}



  if(file.exists("Readme.txt")){

    name <- c(paste0(data.norm,".",project,".matrix"))

    tme <- Sys.Date()
    tme <- format(tme, format="%B %d %Y")


    sp <- data.frame("Type"="---------:" ,"Description"="---------")
    mod <- data.frame("Type"=paste(name, "added the: ") ,"Description"=tme)
    dt <- rbind(sp,mod,sp)


    if(str_detect(name, c("matrix"))==T){
      nr <- nrow(Metadata[[name]])
      nc <- ncol(Metadata[[name]])

      if(str_detect(name, "Raw")){ Assay = "Raw counts"}
      if(str_detect(name, "TPM")){ Assay="TPM normalization"}
      if(str_detect(name, "FKPM")){ Assay="FKPM normalization"}
      if(str_detect(name, "Normalized")){ Assay="Normalized gene expression"}

      ltest <- data.frame("Type"=c("File: " , "Class: " ,"Dimension: ", "Assay: ", "Rownames: ","Colnames: "),"Description" = c(
        paste0(name,".csv"),
        class(Metadata[[name]]),
        paste(nr,"x",nc),
        Assay,
        paste(rownames(Metadata[[name]])[1],"...",rownames(Metadata[[name]])[nrow(Metadata[[name]])]),
        paste(colnames(Metadata[[name]])[2],"...",colnames(Metadata[[name]])[ncol(Metadata[[name]])])
      ))
      }

      dt <- rbind(dt, ltest,sp)
      write.table(dt,"Readme.txt",row.names = F,quote = FALSE, append = T, col.names=FALSE)
      file.show("Readme.txt")
      closeAllConnections()

      }



  return(Metadata) }





    }else{ stop("No meta data object found")} #else meta data object not found











} # function
