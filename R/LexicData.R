#' LexicData function to prepare constant colnames to find in clinical data from lexiue .txt
#'
#' @param list.files.path file path to find lexique of colnames
#' @importFrom data.table fread
#' @import dplyr
#' @importFrom stats na.omit
#' @return a list
#' @export
#'
#' @examples "none"
LexicData <- function(list.files.path){



if(!exists("PatientLexic", mode = "any") | !exists("SamplesLexic", mode = "any") ){

  pos <- 1
  envir = as.environment(pos)



  if(file.exists(paste0(list.files.path$Project.Processes, "/",project,".PatientLexic.txt"))){

    message(paste("Importing PatientLexic.txt from",project, "directory."))
    x <- scan(paste0(list.files.path$Project.Processes, "/",project,".PatientLexic.txt"), what="", sep="\n")%>%strsplit("[[:space:]]+")
    names(x) <- sapply(x, `[[`, 1)
    PL<- lapply(x, `[`, -1)

    assign("PatientLexic", PL, envir = envir)} else {    PL <- scan(paste(list.files.path$Processes,"PatientLexic.txt",sep = "/"), what="", sep="\n")%>%strsplit("[[:space:]]+")
    names(PL) <- sapply(PL, `[[`, 1)
    PL<- lapply(PL, `[`, -1)
    assign("PatientLexic", PL, envir = envir)    }

  if(file.exists(paste0(list.files.path$Project.Processes, "/",project,".SamplesLexic.txt"))){
    message(paste("ImportingSamplesLexic.txt from",project, "directory."))
    x <- scan(paste0(list.files.path$Project.Processes, "/",project,".SamplesLexic.txt"), what="", sep="\n")%>%strsplit("[[:space:]]+")
    names(x) <- sapply(x, `[[`, 1)
    SL<- lapply(x, `[`, -1)

    assign("SamplesLexic", SL, envir = envir)} else {
      SL <- scan(paste(list.files.path$Processes,"SamplesLexic.txt",sep = "/"), what="", sep="\n")%>%strsplit("[[:space:]]+")
      names(SL) <- sapply(SL, `[[`, 1)
      SL<- lapply(SL, `[`, -1)
      assign("SamplesLexic", SL, envir = envir)  }


}

}