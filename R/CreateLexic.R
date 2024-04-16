#' CreateLexic function to prepare constant colnames to find in clinical data from lexiue .txt
#'
#' @param Dataset file path to find lexique of colnames
#' @param type default F. If MetaObject is made of Single.Cell sequencing data, set to T.
#' @importFrom data.table fread
#' @import dplyr
#' @importFrom stats na.omit
#' @return a list
#' @export
#'
#' @examples "none"

CreateLexic <- function(Dataset=NULL, type = c("SamplesLexic","PatientsLexic")){


  if(is.null(Dataset)){stop("Dataset is null.")}
  if(!is.list(Dataset)){stop("Dataset must be a list created from CreateDataset() function.")}
  if(is.null(attributes(Dataset)$File.path) | is.null(attributes(Dataset)$Project)){stop("No attributes found in dataset. Create one with CreateDataset()")}
  if(is.null(attributes(Dataset)$Omics.type)){stop("No Omics.type attributes found in dataset. Firstly add ExpressionMatrix using AddExpressionMatrixRNAseq/Ma/SC functions.")}
  if(!type%in%c("SamplesLexic","PatientsLexic")){stop("Type must be SamplesLexic or PatientsLexic.")}
   list.files.path = attributes(Dataset)$File.path
   project =   attributes(Dataset)$Project
   SC.Lexic=F
   if(attributes(Dataset)$Omics.type=="Single.Cell"){SC.Lexic=T}



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
  if(SC.Lexic==F) {

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
      attr(SL, "Name") = "SC.SamplesLexic"
      return(SL)
      }

  }

  if(SC.Lexic==T) {

    if(file.exists(paste0(list.files.path$Project.Processes, "/",project,"SC.SamplesLexic.txt"))){

      message(paste("Importing SamplesLexic.txt from",project, "directory."))
      x <- scan(paste0(list.files.path$Project.Processes, "/",project,"SC.SamplesLexic.txt"), what="", sep="\n")%>%strsplit("[[:space:]]+")
      names(x) <- sapply(x, `[[`, 1)
      SL<- lapply(x, `[`, -1)
      attr(SL, "Lexic") = "Yes"
      attr(SL, "Name") = "SC.SamplesLexic"
      return(SL)


    } else {

      if(file.exists(paste0(list.files.path$Project.Processes, "/",project,".SamplesLexic.txt"))){
        SL <- scan(paste0(list.files.path$Project.Processes,"/",project,".SamplesLexic.txt"), what="", sep="\n")%>%strsplit("[[:space:]]+")
        names(SL) <- sapply(SL, `[[`, 1)
        SL<- lapply(SL, `[`, -1)

        } else {
      SL <- scan(paste(list.files.path$Processes,"SC.SamplesLexic.txt",sep = "/"), what="", sep="\n")%>%strsplit("[[:space:]]+")
      names(SL) <- sapply(SL, `[[`, 1)
      SL<- lapply(SL, `[`, -1)}
      attr(SL, "Lexic") = "Yes"
      attr(SL, "Name") = "SC.SamplesLexic"
      return(SL)
    }

  }}



}
