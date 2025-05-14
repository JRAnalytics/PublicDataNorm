#' AddClinicFromObject
#'
#' @param Metadata  a Metadata Object to fill
#' @param object an object to add
#' @param setSamplesID.Column a character string : column to fetch SamplesID from Count matrix.
#' @param setPatientID.Column a character string : column to fetch PatientsID from Count matrix.
#' @param setCellsBarcode.Column a character string : column to fetch CellBarcode from Count matrix.
#' @param name name to apply in Metadata object list
#' @param type can be Samples", "Patients" or "Cells", to defiens Data.type attributes
#' @param Export  TRUE or FALSE. If data to be Exported, set T.
#' @param force.replace set as F. T : replace an already object with the same name
#' @import matrixStats
#'
#' @return a Metadata Object
#' @export
#'
#' @examples "none"
#'
AddClinicFromObject  <- function(Metadata,
                           object,
                           setSamplesID.Column = NULL,
                           setPatientID.Column=NULL,
                           setCellsBarcode.Column=NULL,
                           name,
                           type = c("Samples", "Patients", "Cells"),
                           Export = T,
                           force.replace = F){





  if(!is.list(Metadata)) { stop("Metadata object should be a list")}
  if(is.null(object)) { stop("Object is null")}
  if(is.null(name)) { stop("name should be specified, as character")}
  if(!is.character(name)) { stop("name should be specified as character")}
  if(is.null(type)) { stop("Data.type should be specified")}
  if(is.null(Export)) { stop("Export should be specified")}

  if(!type%in%c("Samples", "Patients", "Cells")) { stop("type should be specified from this values c('Samples', 'Patients'or 'Cells') ")}


  if(!all(str_detect(names(Metadata),name)==F)){
    message("An Object with the same name already exist in MetaObject")
    if(force.replace==F){stop("set force.replace==T to subset object.")}}


    l <- length(Metadata)

    if(type=="Samples"){ExpressionMatrixIdColumn = "samplesID"}
    if(type=="Patients"){ExpressionMatrixIdColumn = "patientsID"}
    if(type=="Cells"){ExpressionMatrixIdColumn = "CellsBarcode"}

    if(type=="Samples"){
      if(!is.null(setSamplesID.Column)){

        colnames(object)[ colnames(object)==setSamplesID.Column] = "samplesID"

      if(!is.null(setPatientID.Column)){
        if(setPatientID.Column==setSamplesID.Column){object$patientsID =object$samplesID} else{
        colnames(object)[ colnames(object)==setPatientID.Column] = "patientsID"}}else {stop("setPatientID.Column must be set")}}else{stop("setSamplesID.Column must be set")}}

    if(type=="Patients"){
      if(!is.null(setPatientID.Column)){
        colnames(object)[ colnames(object)==setPatientID.Column] = "patientsID"}else{stop("setPatientID.Column must be set")}
      if(!is.null(setSamplesID.Column)){
        if(setSamplesID.Column==setPatientID.Column){object$samplesID =object$patientsID  } else {colnames(object)[ colnames(object)==setSamplesID.Column] = "samplesID"}}}


    if(type=="Cells"){
      if(!is.null(setCellsBarcode.Column)){colnames(object)[ colnames(object)==setCellsBarcode.Column] = "CellsBarcode"
      }else{stop("setCellsBarcode.Column must be set")}
      if(!is.null(setSamplesID.Column)){   colnames(object)[ colnames(object)==setSamplesID.Column] = "samplesID"}
      if(!is.null(setPatientID.Column)){  colnames(object)[ colnames(object)==setPatientID.Column] = "patientsID"}
      object$CellsBarcode = gsub("[[:punct:]]","-",  object$CellsBarcode)
    }

    if(type != "Patients"){type2 = "Samples annotation"}else{type2 = "Patients clinical data"}
    if(!ExpressionMatrixIdColumn%in%colnames(object)){stop(paste(ExpressionMatrixIdColumn, "is not in",type2, "object colnames."))}




    if(!all(str_detect(names(Metadata),name)==F)){
      if(force.replace==F){stop("set force.replace==T to subset object.")}
      message("Subsetting object.")

      Metadata[[name]] <- object

      t= which(str_detect(names(Metadata),name))

      if(!type%in%c("Samples","Patients", "Cells")){stop("type must be set to Samples, Patients or Cells")}

    if(type == "Samples") {attributes(Metadata)$Data.Type[t] <-  c("SamplesAnnot")}
        if(type == "Patients") {attributes(Metadata)$Data.Type[t] <-  c("Clinic")}
        if(type == "Cells") {attributes(Metadata)$Data.Type[t] <-  c("CellsAnnot")}

        if(Export==T){attributes(Metadata)$Export[t] <- c("Yes") } else {attributes(Metadata)$Export[t] <- c("No") }
      attributes(Metadata)$Cleaned[t] = "No"
      return(Metadata)


    } else {

      Metadata[[l+1]] <- object
      names(Metadata)[l+1] <- name


    if(!type%in%c("Samples","Patients", "Cells")){stop("type must be set to Samples, Patients or Cells")}

    if(l==0) {   if(type == "Samples") {attributes(Metadata)$Data.Type <-  c("SamplesAnnot")}
      if(type == "Patients") {attributes(Metadata)$Data.Type <-  c("Clinic")}
      if(type == "Cells") {attributes(Metadata)$Data.Type <-  c("CellsAnnot")}

      if(Export==T){attributes(Metadata)$Export <- c("Yes") } else {attributes(Metadata)$Export <- c("No") }
      attributes(Metadata)$Cleaned = c( "No")

    } else {  if(type == "Samples") {attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"SamplesAnnot")}
      if(type == "Patients") {attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"Clinic")}
      if(type == "Cells") {attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"CellsAnnot")}

      if(Export==T){attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes") } else {attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"No") }
      attributes(Metadata)$Cleaned = c( attributes(Metadata)$Cleaned,"No")
    }}






  return(Metadata)

}
