#' AddClinic to Meta object
#'
#' @param Metadata Meta object
#' @param path dir path in which the GDC project is saved, or local files are saved
#' @param name if loca=True, names to apply in Metadata object slot
#' @param merge merge loaded data clinic with existing clinial data : full_join by rownames
#' @param name.local.file name file of interest in path directory, could be multiple names. c("a.csv","b.csv")
#' @param mergeBy colname using for merging clinical data.
#' @importFrom utils menu
#' @import purrr
#' @import dplyr
#' @return a data.frame in the Meta Object
#' @export
#'
#' @examples "none"
AddClinic <- function(Metadata, path, merge = c(F,T), mergeBy, name, name.local.file = NULL) {

  ### ecrasement si même nom dans le Meta à faire.


  if(!is.null(Metadata)){

    if(!is.list(Metadata)){
      stop("Metadata should be a list.")}

      l <-length(names(Metadata))

      if(!is.null(name.local.file)){

      filepath <- paste(path,name.local.file,sep="/")

      if(length(filepath)>1){

        if(is.null(mergeBy)){stop("For merging data from multiple loading, mergeBy='colnames' must be specified")}

        clinic <- list()
        count <- 0
        for (i in filepath) {
          if(all(str_detect(i, ".rds|.txt|.csv|.tsv", negate = FALSE)==F) ){stop("No '*.rds' or '*.txt' '.csv' files in set directory. \n change path or add file")}
          count <- count+1

          if(str_detect(i, ".rds", negate = FALSE)){
            clinic[[count]] <- readRDS(i) }

          if(str_detect(i, ".txt|.csv|.tsv", negate = FALSE)){

            clinic[[count]] <- suppressWarnings(as.data.frame(data.table::fread(i, na.strings = "")))
            rownames(clinic[[count]]) <-    clinic[[count]][,1]

          }}

          dt <- clinic %>% purrr::reduce(full_join, by=mergeBy)
          rownames(dt) <- dt[,1]


        }
      if(length(filepath)==1){

        if(all(str_detect(filepath, ".rds|.txt|.csv|.tsv", negate = FALSE)==F) ){stop("No '*.rds' or '*.txt' '.csv' files in set directory. \n change path or add file")}


        if(str_detect(filepath, ".rds", negate = FALSE)){
          dt <- readRDS(filepath) }

        if(str_detect(filepath, ".txt|.csv|.tsv", negate = FALSE)){

          dt <- suppressWarnings(as.data.frame(data.table::fread(filepath, na.strings = "")))
          rownames(dt) <-    dt[,1]

        }



      }} else {

          LF <- list.files(path)
          LF <- LF[str_detect(LF,"clinic")]
          if(length(LF)>1){stop("Ther is more than one files with 'clinic' in its name. Switch to 'name.local.file' and 'mergeBy' ")}

      if(all(str_detect(LF, ".rds|.txt|.csv|.tsv", negate = FALSE)==F)  ){stop("#55 No '*.rds' or '*.txt' '.csv' files in set directory. \n change path or add file")}


      if(str_detect(LF, ".rds", negate = FALSE)){Metadata[[l+1]] <- readRDS(LF) }

      if(str_detect(LF, ".txt|.csv|.tsv", negate = FALSE)){

        dt <- suppressWarnings(as.data.frame(data.table::fread(LF, na.strings = "")))
        rownames(dt) <- dt[,1]
        }}


        if (merge == F){


            Metadata[[l+1]] <- dt

            names(Metadata)[l+1] <- name
        } else {

          if(is.null(str_detect(names(Metadata),"clinic"))){stop("No clinical data found in Meta Object.\n
                                                                 For merging 1 newly loaded clinical data to an already loaded file in Meta Object, firslty load one.\n
                                                                 You can list a list of caracter for multiple clinical data to load at once and will be full_join.")}

          if(is.null(mergeBy)){stop("For merging data, mergeBy='colnames' must be specified")}
          clinic <- list(Metadata[which(str_detect(names(Metadata),"clinic"))][[1]],dt)

          clinic <- clinic %>% purrr::reduce(full_join, by=mergeBy)

          Metadata[which(str_detect(names(Metadata),"clinic"))][[1]] <- clinic
        }

      return(Metadata)


    }else{ stop("No meta data object found")

  } } # function
