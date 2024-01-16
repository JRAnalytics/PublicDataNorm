#' AddClinicFromFile to Meta object
#'
#' @param Metadata Meta object
#' @param Raw.file.path dir path in which the GDC project is saved, or local files are saved
#' @param name if loca=True, names to apply in Metadata object slot
#' @param type c("Samples", "Patients") To specify if clinical data refers to samples or patients data. Importat if mutliple samples for one patient.
#' @param mergeToClinic merge loaded data clinic with existing clincial data : full_join by rownames. False if you load more than 1 file.
#' @param name.local.file name file of interest in path directory, could be multiple names. c("a.csv","b.csv")
#' @param mergeBy colname using for merging clinical data.
#' @param Export  TRUE or FALSE. If data to be Exported, set T.
#' @param join c("left_join", "full_join")
#' @param force.replace set as F. T : replace an already object with the same name
#' @importFrom utils menu
#' @import purrr
#' @import dplyr
#' @return a data.frame in the Meta Object
#' @export
#'
#' @examples "none"
AddClinicFromFile <- function(Metadata,
                              Raw.file.path,
                              name.local.file = NULL,
                              name,
                              type = c("Samples", "Patients", "Cells"),
                              force.replace=F,
                              Export = F,
                              mergeToClinic = F,
                              mergeBy,
                              join = c("left_join", "full_join")) {

  ### ecrasement si même nom dans le Meta à faire.


  path = Raw.file.path

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

          }}

        if(join=="full_join"){
          dt <- clinic %>% purrr::reduce(full_join, by=mergeBy)
        }

        if(join=="left_join"){
          dt <- clinic %>% purrr::reduce(left_join, by=mergeBy)
        }




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


        if(str_detect(LF, ".rds", negate = FALSE)){
          dt <- readRDS(LF)}

        if(str_detect(LF, ".txt|.csv|.tsv", negate = FALSE)){

          dt <- suppressWarnings(as.data.frame(data.table::fread(LF, na.strings = "")))
          rownames(dt) <- dt[,1]
        }}


    if (mergeToClinic == F){


      if(!all(str_detect(names(Metadata),name)==F)){
        message("An Object with the same name already exist in MetaObject")
        if(force.replace==F){stop("set force.replace==T to subset object.")}
        message("Subsetting object.")
        Metadata[[name]] <- dt

        tt = which(str_detect(names(Metadata), name))

        if(!type%in%c("Samples","Patients")){stop("type must be set to Samples or Patients")}

         if(type == "Samples") {attributes(Metadata)$Data.Type[tt] <-  "SamplesAnnot"}
          if(type == "Patients") {attributes(Metadata)$Data.Type[tt] <- "Clinic"}
          if(Export==T){attributes(Metadata)$Export[tt] <- "Yes" } else {
            attributes(Metadata)$Export[tt] <- "No" }


      } else { Metadata[[l+1]] <- dt
      names(Metadata)[l+1] <- name



      if(!type%in%c("Samples","Patients")){stop("type must be set to Samples or Patients")}

      if(l==0) {   if(type == "Samples") {attributes(Metadata)$Data.Type <-  c("SamplesAnnot")}
        if(type == "Patients") {attributes(Metadata)$Data.Type <-  c("Clinic")}
        if(type == "Cells") {attributes(Metadata)$Data.Type <-  "CellsAnnot"}
        if(Export==T){attributes(Metadata)$Export <- c("Yes") } else {attributes(Metadata)$Export <- c("No") }


      } else {  if(type == "Samples") {attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"SamplesAnnot")}
        if(type == "Patients") {attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"Clinic")}
        if(type == "Cells") {attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"CellsAnnot")}
        if(Export==T){attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes") } else {attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"No") }

      }}
    } else {

      if(is.null(str_detect(attributes(Metadata)$Data.Type,"Clinic"))){stop("No clinical data found in Meta Object.\n
                                                                 For merging 1 newly loaded clinical data to an already loaded file in Meta Object, firslty load one.\n
                                                                 You can list a list of caracter for multiple clinical data to load at once and will be full_join.")}

      if(is.null(mergeBy)){stop("For merging data, mergeBy='colnames' must be specified")}

      NB <- which(str_detect(string =attributes(Metadata)$Data.Type ,"Clinic") & attributes(Metadata)$Export=="No")

      clinic <- list(Metadata[[NB]], dt)

      if(join=="full_join"){
        dt <- clinic %>% purrr::reduce(full_join, by=mergeBy)
      }

      if(join=="left_join"){
        dt <- clinic %>% purrr::reduce(left_join, by=mergeBy)
      }

      Metadata[[NB]] <- clinic


    }

    return(Metadata)


  }else{ stop("No meta data object found")

  } } # function
