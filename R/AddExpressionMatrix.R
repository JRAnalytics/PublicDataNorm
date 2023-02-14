#' AddExpressionMatrix to Meta object
#'
#' @param Metadata Meta object
#' @param query query file from TCGAimportExpression function
#' @param data.norm a character "Raw" ,"TPM", "FKPM"
#' @param path dir path in which the GDC project is saved, or local files are saved
#' @param local T of F. If F, use query form. If F, add expression patrice from local file
#' @param name.local.file if loca=True, names to apply in Meatadata object slot
#' @importFrom utils menu
#' @import data.table
#' @return a data.frame in the Meta Object
#' @export
#'
#' @examples "none"
AddExpressionMatrix <- function(Metadata, local = c(T, F) , query, data.norm, path, name.local.file) {
  if(!is.null(Metadata)){
    if(!is.list(Metadata)){
      stop("Metadata should be a list.")}

  setwd(path)
  if(local== T){




     message("Local import")

    l <-length(names(Metadata))
    lf <- list.files(path)

    if(length(lf)>1){print(c(message("There is more than one Matrix files in Dir :"),lf))}

    if(all(str_detect(lf, ".rds|.txt|.csv|.tsv", negate = FALSE)==F)){stop("No '*.rds' or '*.txt' or '*.csv' files in set directory. \n change path or add file")}

    if(length(lf[str_detect(lf, "matrix")])>1) {

      for (i in lf[str_detect(lf, "matrix")]) {
        l <- length(Metadata)
        message(paste("Loading", i, "file"))

        if(str_detect(i, ".rds", negate = FALSE)){Metadata[[1]] <- readRDS(i)} else {
          if(str_detect(i, ".txt", negate = FALSE)){

            dt <- suppressWarnings(as.data.frame(data.table::fread(i)))
            rownames(dt) <- dt[,1]
            Metadata[[l+1]] <- dt}  else {
              if(str_detect(i, ".csv", negate = FALSE)){
                dt <- suppressWarnings(as.data.frame(data.table::fread(i)))
                rownames(dt) <- dt[,1]
                Metadata[[l+1]] <- dt} else {

                  if(str_detect(i, ".tsv", negate = FALSE)){
                    dt <- suppressWarnings(as.data.frame(data.table::fread(i)))
                    rownames(dt) <- dt[,1]
                    Metadata[[l+1]] <- dt}}
            } }

        names(Metadata)[l+1] <- paste0(name.local.file,".matrix.",which(lf[str_detect(lf, "matrix")]%in%i))

        } } else {

    if(length(Metadata)>1) {


      lf <- lf[str_detect(lf, "matrix")]
      message(paste("Loading", lf, "file"))

    if(str_detect(lf, ".rds", negate = FALSE)){Metadata[[l+1]] <- readRDS(lf)} else {

      if(str_detect(lf, ".txt", negate = FALSE)){
        dt <- suppressWarnings(as.data.frame(data.table::fread(lf)))
        rownames(dt) <- dt[,1]
        dt <- dt[,colnames(Metadata[[1]])]
        Metadata[[l+1]] <- dt}  else {

          if(str_detect(lf, ".csv", negate = FALSE)){
            dt <- suppressWarnings(as.data.frame(data.table::fread(lf)))
            rownames(dt) <- dt[,1]
            dt <- dt[,colnames(Metadata[[1]])]
            Metadata[[l+1]] <- dt} else {

              if(str_detect(lf, ".tsv", negate = FALSE)){
              dt <- suppressWarnings(as.data.frame(data.table::fread(lf)))
              rownames(dt) <- dt[,1]
              dt <- dt[,colnames(Metadata[[1]])]
              Metadata[[l+1]] <- dt}
              }
          }

      }

    names(Metadata)[l+1] <- paste0(name.local.file,".matrix")}

     else {

       lf <- lf[str_detect(lf, "matrix")]
       message(paste("Loading", lf, "file"))
      if(str_detect(lf, ".rds", negate = FALSE)){Metadata[[1]] <- readRDS(lf)} else {
        if(str_detect(lf, ".txt", negate = FALSE)){

        dt <- suppressWarnings(as.data.frame(data.table::fread(lf)))
        rownames(dt) <- dt[,1]
        Metadata[[1]] <- dt}  else {
          if(str_detect(lf, ".csv", negate = FALSE)){
            dt <- suppressWarnings(as.data.frame(data.table::fread(lf)))
            rownames(dt) <- dt[,1]
            Metadata[[1]] <- dt}else {

              if(str_detect(lf, ".tsv", negate = FALSE)){
                dt <- suppressWarnings(as.data.frame(data.table::fread(lf)))
                rownames(dt) <- dt[,1]
                Metadata[[1]] <- dt}}
        }
        }

      names(Metadata)[1] <- paste0(name.local.file,".matrix")




      }}


###marche pas le readme!

    if(file.exists("Readme.txt")){


      name <- paste0(name.local.file,".matrix")

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


    return(Metadata)


  } else { # if local ==T
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
    genome = ifelse(query$legacy,"hg19","hg38"),
    experimental.strategy = unique(query$results[[1]]$experimental_strategy))

  y <- y[,colnames(Metadata[[1]])]

  if(!all(colnames(Metadata[[1]])==colnames(y))){stop("Samples are not the sames accross data matrix expression")}


  l <-length(names(Metadata))
  Metadata[[l+1]]<- y
  names(Metadata)[l+1] <- c(paste0(data.norm,".",project,".matrix"))



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
