#' AddClinicFromFile to Meta object
#'
#' @param Metadata Meta object
#' @param name if loca=True, names to apply in Metadata object slot
#' @param type c("Samples", "Patients") To specify if clinical data refers to samples or patients data. Importat if mutliple samples for one patient.
#' @param mergeToClinic default NULL. name of loaded clinical data in Metadata. Merge the newly added clinical data to a already loaded data clinic with existing clincial data : full_join by rownames.
#' @param ClinicFile name file of interest in path directory, could be multiple names. c("a.csv","b.csv")
#' @param ExpressionMatrixIdColumn column to fetch colnames from Count matrix.
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
                              ClinicFile = NULL,
                              ExpressionMatrixIdColumn=NULL,
                              name= NULL,
                              type = c("Samples", "Patients", "Cells"),
                              force.replace=F,
                              Export = F,
                              mergeToClinic = NULL,
                              mergeBy= NULL,
                              join = c("left_join", "full_join")) {



  if(is.null(Metadata)){stop("No meta data object found")}
  if(!is.list(Metadata)){stop("Metadata should be a list.")}
  if(is.null(ClinicFile)){stop("No set ClinicFile information")}
  if(!inherits(ClinicFile,"character")){stop("ClinicFile is not a character string.") }
  if(is.null(ExpressionMatrixIdColumn)){stop("ExpressionMatrixIdColumn must be specifierdt (i.e character string : colname of clinic to add).")}
  if(length(ExpressionMatrixIdColumn)>1){stop("ExpressionMatrixIdColumn must be of length equal to one.")}

  l <-length(names(Metadata))
  filepath <- paste(attributes(Metadata)$File.path$Project.RawData,ClinicFile,sep="/")


  if(length(filepath)>1){

        if(is.null(mergeBy)){stop("For merging data from multiple loading, mergeBy='colnames' must be specified")}

        clinic <- list()
        count <- 0
        for (i in filepath) {
          if(all(str_detect(i, ".rds|.txt|.csv|.tsv", negate = FALSE)==F) ){stop("No '*.rds' or '*.txt' '.csv' files in set directory. \n change path or add file")}
          count <- count+1

          if(str_detect(i, ".rds", negate = FALSE)){clinic[[count]] <- readRDS(i) }

          if(str_detect(i, ".txt|.csv|.tsv", negate = FALSE)){clinic[[count]] <- suppressWarnings(as.data.frame(data.table::fread(i, na.strings = "")))}}

        if(join=="full_join"){dt <- clinic %>% purrr::reduce(full_join, by=mergeBy)}

        if(join=="left_join"){dt <- clinic %>% purrr::reduce(left_join, by=mergeBy)}}





      if(length(filepath)==1){

        if(all(str_detect(filepath, ".rds|.txt|.csv|.tsv", negate = FALSE)==F) ){stop("No '*.rds' or '*.txt' '.csv' files in set directory. \n change path or add file")}


        if(str_detect(filepath, ".rds", negate = FALSE)){dt <- readRDS(filepath) }

        if(str_detect(filepath, ".txt|.csv|.tsv", negate = FALSE)){
          dt <- suppressWarnings(as.data.frame(data.table::fread(filepath, na.strings = "")))
          rownames(dt) <-    dt[,1]}
        }


    if (is.null(mergeToClinic)){


      if(!all(str_detect(names(Metadata),name)==F)){
        message("An Object with the same name already exist in MetaObject")
        if(force.replace==F){stop("set force.replace==T to subset object.")}
        message("Subsetting object.")
        Metadata[[name]] <- dt

        tt = which(str_detect(names(Metadata), name))

        if(!type%in%c("Samples","Patients","Cells")){stop("type must be set to Samples or Patients")}

         if(type == "Samples") {attributes(Metadata)$Data.Type[tt] <-  "SamplesAnnot"}
          if(type == "Patients") {attributes(Metadata)$Data.Type[tt] <- "Clinic"}
        if(type == "Cells") {attributes(Metadata)$Data.Type[tt] <- "CellsAnnot"}
          if(Export==T){attributes(Metadata)$Export[tt] <- "Yes" } else {
            attributes(Metadata)$Export[tt] <- "No" }


      } else { Metadata[[l+1]] <- dt
      names(Metadata)[l+1] <- name



      if(!type%in%c("Samples","Patients", "Cells")){stop("type must be set to Samples or Patients")}

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

      if(!mergeToClinic%in%names(Metadata)){stop(paste("No",mergeToClinic,"found in Meta Object.\n
                                                                 For merging 1 newly loaded clinical data to an already loaded file in Meta Object, firslty load one.\n
                                                                 You can list a list of caracter for multiple clinical data to load at once and will be full_join."))}

      if(is.null(mergeBy)){stop("For merging data, mergeBy='colnames' must be specified")}


      dt <- list(Metadata[[mergeToClinic]], dt)

      if(join=="full_join"){
        dt <- dt %>% purrr::reduce(full_join, by=mergeBy)
      }

      if(join=="left_join"){
        dt <- dt %>% purrr::reduce(left_join, by=mergeBy)
      }

      Metadata[[mergeToClinic]] <- dt


    }

if(!ExpressionMatrixIdColumn%in%colnames(dt)){stop(paste(ExpressionMatrixIdColumn, "colnames is not in the added clinical data."))}


if(type!="Cells"){
  zz <- which(attributes(Metadata)$Data.Type=="Count")[1]
  samples <- colnames(Metadata[[zz]])
  message("Found Samples :")
  print(summary(samples%in%dt[,ExpressionMatrixIdColumn]))
  message("Unfound Samples :")
  print(samples[!samples%in%dt[,ExpressionMatrixIdColumn]])}


    return(Metadata)



  } # function
