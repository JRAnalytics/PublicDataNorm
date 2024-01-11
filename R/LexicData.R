#' LexicData function to prepare constant colnames to find in clinical data from lexiue .txt
#'
#' @param list.files.path file path to find lexique of colnames
#' @param replaceLexic replace Samples and Patient Lexis if already loaded
#' @param SC.Lexic default F. If MetaObject is made of Single.Cell sequencing data, set to T.
#' @importFrom data.table fread
#' @import dplyr
#' @importFrom stats na.omit
#' @return a list
#' @export
#'
#' @examples "none"

LexicData <- function(list.files.path, SC.Lexic = F, replaceLexic = F){



if(c(exists("PatientLexic", mode = "any") | exists("SamplesLexic", mode = "any")) & replaceLexic==T |
   c(!exists("PatientLexic", mode = "any") | !exists("SamplesLexic", mode = "any")) ){



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


  if(SC.Lexic==F) {

    if(file.exists(paste0(list.files.path$Project.Processes, "/",project,".SamplesLexic.txt"))){

      message(paste("Importing SamplesLexic.txt from",project, "directory."))
      x <- scan(paste0(list.files.path$Project.Processes, "/",project,".SamplesLexic.txt"), what="", sep="\n")%>%strsplit("[[:space:]]+")
      names(x) <- sapply(x, `[[`, 1)
      SL<- lapply(x, `[`, -1)

      assign("SamplesLexic", SL, envir = envir)


    } else {

      SL <- scan(paste(list.files.path$Processes,"SamplesLexic.txt",sep = "/"), what="", sep="\n")%>%strsplit("[[:space:]]+")
      names(SL) <- sapply(SL, `[[`, 1)
      SL<- lapply(SL, `[`, -1)

      assign("SamplesLexic", SL, envir = envir)
      }

  }

  if(SC.Lexic==T) {

    if(file.exists(paste0(list.files.path$Project.Processes, "/",project,"SC.SamplesLexic.txt"))){

      message(paste("Importing SamplesLexic.txt from",project, "directory."))
      x <- scan(paste0(list.files.path$Project.Processes, "/",project,"SC.SamplesLexic.txt"), what="", sep="\n")%>%strsplit("[[:space:]]+")
      names(x) <- sapply(x, `[[`, 1)
      SL<- lapply(x, `[`, -1)

      assign("SamplesLexic", SL, envir = envir)


    } else {

      if(file.exists(paste0(list.files.path$Project.Processes, "/",project,".SamplesLexic.txt"))){
        SL <- scan(paste(list.files.path$Processes,"SamplesLexic.txt",sep = "/"), what="", sep="\n")%>%strsplit("[[:space:]]+")
        names(SL) <- sapply(SL, `[[`, 1)
        SL<- lapply(SL, `[`, -1)

        } else {
      SL <- scan(paste(list.files.path$Processes,"SC.SamplesLexic.txt",sep = "/"), what="", sep="\n")%>%strsplit("[[:space:]]+")
      names(SL) <- sapply(SL, `[[`, 1)
      SL<- lapply(SL, `[`, -1)}

      assign("SamplesLexic", SL, envir = envir)
    }

  }





} else {  stop("Patientlexic and SamplesLexic already loaded. Set replaceLexic=T to replace lexics.")}

}
