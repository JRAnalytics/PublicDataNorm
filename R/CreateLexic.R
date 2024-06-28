#' CreateLexic function to prepare constant colnames to find in clinical data from lexiue .txt
#'
#' @param Metadata file path to find lexique of colnames
#' @param type c("SamplesLexic","PatientsLexic","CellsLexic")
#' @importFrom data.table fread
#' @import dplyr
#' @importFrom stats na.omit
#' @return a list
#' @export
#'
#' @examples "none"

CreateLexic <- function(Metadata=NULL, type = c("SamplesLexic","PatientsLexic","CellsLexic")){


  if(is.null(Metadata)){stop("Metadata is null.")}
  if(!is.list(Metadata)){stop("Metadata must be a list created from CreateMetadata() function.")}
  if(is.null(attributes(Metadata)$File.path) | is.null(attributes(Metadata)$Project)){stop("No attributes found in Metadata. Create one with CreateMetadata()")}
  if(is.null(attributes(Metadata)$Omics.type)){stop("No Omics.type attributes found in Metadata. Firstly add ExpressionMatrix using AddExpressionMatrixRNAseq/Ma/SC functions.")}
  if(!type%in%c("SamplesLexic","PatientsLexic", "CellsLexic")){stop("Type must be SamplesLexic or PatientsLexic.")}
   list.files.path = attributes(Metadata)$File.path
   project =  attributes(Metadata)$Project




  if(type=="PatientsLexic"){
  if(file.exists(paste0(list.files.path$Project.Processes, "/",project,".PatientLexic.txt"))){

    message(paste("Importing PatientLexic.txt from",project, "directory."))
    x <- scan(paste0(list.files.path$Project.Processes, "/",project,".PatientLexic.txt"), what="", sep="\n")%>%strsplit("[[:space:]]+")
    names(x) <- sapply(x, `[[`, 1)
    PL<- lapply(x, `[`, -1)

    attr(PL, "Lexic") = "Yes"
    attr(PL, "Name") = "PatientsLexic"

    return(PL)} else {    PL <- scan(paste(list.files.path$Processes,"PatientLexic.txt",sep = "/"), what="", sep="\n")%>%strsplit("[[:space:]]+")
    names(PL) <- sapply(PL, `[[`, 1)
    PL<- lapply(PL, `[`, -1)

    attr(PL, "Lexic") = "Yes"
    attr(PL, "Name") = "PatientsLexic"

     return(PL)  }}

  if(type=="SamplesLexic"){


    if(file.exists(paste0(list.files.path$Project.Processes, "/",project,".SamplesLexic.txt"))){

      message(paste("Importing SamplesLexic.txt from",project, "directory."))
      x <- scan(paste0(list.files.path$Project.Processes, "/",project,".SamplesLexic.txt"), what="", sep="\n")%>%strsplit("[[:space:]]+")
      names(x) <- sapply(x, `[[`, 1)
      SL<- lapply(x, `[`, -1)

      attr(SL, "Lexic") = "Yes"
      attr(SL, "Name") = "SamplesLexic"

      return(SL)


    } else {

      SL <- scan(paste(list.files.path$Processes,"SamplesLexic.txt",sep = "/"), what="", sep="\n")%>%strsplit("[[:space:]]+")
      names(SL) <- sapply(SL, `[[`, 1)
      SL<- lapply(SL, `[`, -1)
      attr(SL, "Lexic") = "Yes"
      attr(SL, "Name") = "SamplesLexic"
      return(SL)


  }
}


   if(type=="CellsLexic"){


       if(file.exists(paste0(list.files.path$Project.Processes, "/",project,"CellsLexic.txt"))){

         message(paste("Importing SamplesLexic.txt from",project, "directory."))
         x <- scan(paste0(list.files.path$Project.Processes, "/",project,"CellsLexic.txt"), what="", sep="\n")%>%strsplit("[[:space:]]+")
         names(x) <- sapply(x, `[[`, 1)
         SL<- lapply(x, `[`, -1)
         attr(SL, "Lexic") = "Yes"
         attr(SL, "Name") = "CellsLexic"
         return(SL)


       } else {

         if(file.exists(paste0(list.files.path$Project.Processes, "/",project,".CellsLexic.txt"))){
           SL <- scan(paste0(list.files.path$Project.Processes,"/",project,".CellsLexic.txt"), what="", sep="\n")%>%strsplit("[[:space:]]+")
           names(SL) <- sapply(SL, `[[`, 1)
           SL<- lapply(SL, `[`, -1)

         } else {
           SL <- scan(paste(list.files.path$Processes,"CellsLexic.txt",sep = "/"), what="", sep="\n")%>%strsplit("[[:space:]]+")
           names(SL) <- sapply(SL, `[[`, 1)
           SL<- lapply(SL, `[`, -1)}
         attr(SL, "Lexic") = "Yes"
         attr(SL, "Name") = "CellsLexic"
         return(SL)
      }}



}
