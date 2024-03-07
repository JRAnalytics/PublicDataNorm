#' AddClinicFromObject
#'
#' @param Metadata  a Metadata Object to fill
#' @param object an object to add
#' @param ExpressionMatrixIdColumn column to fetch colnames from Count matrix.
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
                           name,
                           ExpressionMatrixIdColumn = NULL,
                           type = c("Samples", "Patients", "Cells"),
                           Export = T,
                           force.replace = F){



  if(!is.list(Metadata)) { stop("Metadata object should be a list")}
  if(is.null(object)) { stop("Object is null")}
  if(is.null(name)) { stop("name should be specified, as character")}
  if(!is.character(name)) { stop("name should be specified as character")}
  if(is.null(type)) { stop("Data.type should be specified")}
  if(is.null(Export)) { stop("Export should be specified")}
  #if(is.null(ExpressionMatrixIdColumn)){stop("ExpressionMatrixIdColumn must be specifierdt (i.e character string : colname of clinic to add).")}
  if(!type%in%c("Samples", "Patients", "Cells")) { stop("type should be specified from this values c('Samples', 'Patients'or 'Cells') ")}



  if(!all(str_detect(names(Metadata),name)==F)){
    message("An Object with the same name already exist in MetaObject")
    if(force.replace==F){stop("set force.replace==T to subset object.")}}


    l <- length(Metadata)

    if(type=="Samples"){
    object =as.data.frame(object)
    zz <- which(attributes(Metadata)$Data.Type=="Count")[1]
    samples <- colnames(Metadata[[zz]])
    message("Found Samples :")
    print(summary(samples%in%object[,ExpressionMatrixIdColumn]))
    message("Unfound Samples :")
    print(samples[!samples%in%object[,ExpressionMatrixIdColumn]])}

    # if(SamplesFilter==T){
    #   message("Selecting only Samples present in both Count and Clinical data.")
    #
    #   # Diviser chaque élément de la colonne en un vecteur de sous-chaînes
    #   substrings <- strsplit(as.data.frame(object)[,colT], ";")
    #   # Vérifier si chaque valeur du vecteur est présente dans chaque vecteur de sous-chaînes
    #   est_present <- sapply(substrings, function(x) any(colnames(Metadata[[zz]]) %in% x))
    #
    #   # Récupérer les index de position correspondants
    #   indices <- which(est_present)
    #   object <- object[indices,]
    #
    #
    # }


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

      return(Metadata)


    } else {

      Metadata[[l+1]] <- object
      names(Metadata)[l+1] <- name


    if(!type%in%c("Samples","Patients", "Cells")){stop("type must be set to Samples, Patients or Cells")}

    if(l==0) {   if(type == "Samples") {attributes(Metadata)$Data.Type <-  c("SamplesAnnot")}
      if(type == "Patients") {attributes(Metadata)$Data.Type <-  c("Clinic")}
      if(type == "Cells") {attributes(Metadata)$Data.Type <-  c("CellsAnnot")}

      if(Export==T){attributes(Metadata)$Export <- c("Yes") } else {attributes(Metadata)$Export <- c("No") }


    } else {  if(type == "Samples") {attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"SamplesAnnot")}
      if(type == "Patients") {attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"Clinic")}
      if(type == "Cells") {attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"CellsAnnot")}

      if(Export==T){attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes") } else {attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"No") }

    }}






  return(Metadata)

}
